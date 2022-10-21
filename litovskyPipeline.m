clear all
selpath = uigetdir;
addpath(selpath);
addpath('Homer3')
files = dir([selpath '\*.snirf'])
% outpath = 'C:\Users\micza\Documents\SMARTLab\fNIRS\SPIN2021_Converted Data\snirf\litovsky pipeline\processedSnirf';

sheets = sheetnames('excludedChannels.xlsx'); % Read the sheetnames with excel file containing excluded channels based on qtnirs

for i = 1:size(files,1)
    snirf0 = SnirfLoad(files(i).name);

    %% Step 1: Remove step-like noise (applied on light intensity data)
    snirf1 = SnirfClass();
    snirf1 = snirf0;
    intensity = snirf1.data.dataTimeSeries(2:end,:); %Remove the first sample -- it's a spike resulting from device turning on
    I = removeStep(intensity);
    snirf1.data.dataTimeSeries = [I(1,:); I];

    %% Step 2: Exclude 'poor' channels
    snirf2 = snirf1;
    underscoreIdx = strfind(files(i).name,'_');
    snirfIdx = strfind(files(i).name,'.snirf');
    pid{i} = files(i).name(underscoreIdx+1:snirfIdx-1);

    %% Step 3: Convert light intensity to optical density
    snirf3 = snirf2;
    dod = hmrR_Intensity2OD(snirf1);
    snirf3.data = dod;


    %% Step 4: Perform motion correction using the wavelet decomposition
    % method
    data_dod = SnirfClass();
    data_dod = snirf3;
    mlActMan = [];
    mlActAuto = [];
    iqr=1.5;
    turnon=1;

    data_dod = hmrR_MotionCorrectWavelet(data_dod, mlActMan, mlActAuto, iqr, turnon);

    %% Step 5: Convert OD to Concentration
    ppf = 1;
    dc = hmrR_OD2Conc( data_dod.data, data_dod.probe, [ppf ppf] );
    SDpairs = dc.GetMeasListSrcDetPairs; %Get the source-detector pairs
    for chan = 1:size(SDpairs,1)
        chLabels{i}{chan} = ['S',num2str(SDpairs(chan,1)),'-D',num2str(SDpairs(chan,2))];
    end

    %% Step 6: Apply a bandpass filter (0.01-1.5Hz)
    hpf = 0.01;
    lpf = 1.5;
    dc2 = hmrR_BandpassFilt(dc, hpf, lpf);

    %% Step 7: Subtract the short-channel component
    % Get the SSCs
    measList = dc2.GetMeasList;
    srcList{i} = measList(1:3:end,1);

    % The indices for the SSCs (HbO, HbR, HbT)
    rSSCidx = find(measList(:,1)==2);
%     rssc(i) = rSSCidx(1); % Based on testing, this doesn't change
    lSSCidx = find(measList(:,1)==7);
%     lssc(i) = lSSCidx(1); % This changes though.. but we account for it
%     in the boolean statement (==) used above

    % The timeseries data for the SSCs (HbO, HbR, HbT)
    rightSSC = dc2.dataTimeSeries(:,rSSCidx);
    leftSSC = dc2.dataTimeSeries(:,lSSCidx);

    % Get the indices for left/right optodes
    rOptIdx = 1:size(measList,1)/2;
    lOptIdx = size(measList,1)/2+1:size(measList,1);

    %% For HbO, HbR, HbT
    cnt = 1;
    cnt2 = 1;
    for ii = 1:size(measList,1)
        % Identify the predictor variable (X) based on hemisphere
        if ii <= rOptIdx(end) %right hemisphere
            X = rightSSC(:,cnt);
        else %left hemisphere
            X = leftSSC(:,cnt);
        end
        cnt = cnt+1;
        if cnt > 3
            cnt = 1;
        end

        % Get the response variable (y)
        y = dc2.dataTimeSeries(:,ii);

        % Fit GLM to the data HbO ~ SSC and get the beta value
        b = glmfit(X,y,[]);

        % Subtract product of beta and SSC from the response variable
        dcpost(:,cnt2) = y - b(2).*X;
        cnt2 = cnt2+1;
    end
    dc3 = dc2;
    dc3.SetDataTimeSeries(dcpost)
    clear dcpost

    %% Step 8: Apply a third-order butterworth filter (0.01-0.09 Hz)
    hpf = 0.01;
    lpf = 0.09;
    dc4 = hmrR_BandpassFilt(dc3, hpf, lpf);

    %% Step 9: Calculate the block-average for HbO and HbR
    stim = data_dod.stim; %The stim marker
    HbOBlockAvg{i} = blockAvg(dc4, stim, SDpairs); %The block averages

end

%% Remove excluded channels and SSCs for visualizations and stats stuff
HbO_exCh = HbOBlockAvg;  % Copy the HbOBlockAvg into a new variable
HbO_exCh2 = HbOBlockAvg; % For Hannah
HbO_exCh3 = HbOBlockAvg; % For testing modifications to code

% Source indices corresponding to ROIs
RTsrcIdx = [1,3];
RFsrcIdx = [4,5];
LTsrcIdx = [6,8];
LFsrcIdx = [9,10];

RejPid = {}; %Preallocate memory for a list of rejected participants
RejIdx = [];
pidList = {};
pidAcceptIdx = [];
count = 1;

% iii is the index for participants
for iii = 1:size(HbO_exCh,2)
    %     figure %For plotting waveforms for individual participants
   %% Set the excluded channels to NaN values
    pidIdx = find(strcmp(pid{iii},sheets));
    excludedCh = readtable('excludedChannels.xlsx',sheet = sheets(pidIdx));

    % Loop through the list of excluded channels
    for exCh = 1:size(excludedCh,1)
        ch = ['S', num2str(excludedCh{exCh,"Var1"}),'-D', num2str(excludedCh{exCh,"Var2"})]; %Transform numbers into ch labels
        chIdx = find(strcmp(ch,chLabels{iii})); %Find the index of excluded channels

        % Loop through all conditions.
        % Excluded channels will become NaN values
        for cndtns = 1:size(HbO_exCh{iii},2) 
            HbO_exCh{iii}{cndtns}(:,chIdx) = nan;
            HbO_exCh2{iii}{cndtns}(:,chIdx) = nan; %For hannah
%             HbO_exCh3{iii}{cndtns}(:,chIdx) = nan; %For testing
        end
    end

   %% Calculate the mean waveform for every 4 columns (each cluster)
    % jjj is the index for conditions
    for jjj = 1:size(HbOBlockAvg{iii},2)
        %HbOROI{iii,jjj}
        sig = HbO_exCh{iii}{jjj}; %Get the signal
        sig_temp = [sig(:,ismember(srcList{iii},RTsrcIdx)) sig(:,ismember(srcList{iii},RFsrcIdx))...
            sig(:,ismember(srcList{iii},LTsrcIdx)) sig(:,ismember(srcList{iii},LFsrcIdx))];
%         sig3 = HbO_exCh3{iii}{jjj}; %For testing
       
%         % Remove the SSCs
%         leftSC = find(srcList{iii}==7); %Remove the left SSC first -- it's located near end of matrix
%         sig(:,leftSC) = [];
%         rightSC = find(srcList{iii}==2); %Remove the right SSC second -- it's located near beginning of matrix
%         sig(:,rightSC) = [];
%         sig(:,3) = [];
%         sig(:,11) = [];

        %% Check to see if we need to exclude the participant
        nanRejThresh = 3; %Reject the participant if there are 3 NaNs in a cluster
        if any(sum(isnan(reshape(sig_temp(1,:)',[4,4]))) >= nanRejThresh) %Check to see if we need to exclude participant from analysis
            RejPid = [RejPid; pid{iii}];
            RejIdx = [RejIdx; iii];
            break
        else
            % Cluster into ROIs and calculate the mean waveform for the cluster
            RT = mean(sig(:,ismember(srcList{iii},RTsrcIdx)),2,'omitnan');
            RF = mean(sig(:,ismember(srcList{iii},RFsrcIdx)),2,'omitnan');
            LT = mean(sig(:,ismember(srcList{iii},LTsrcIdx)),2,'omitnan');
            LF = mean(sig(:,ismember(srcList{iii},LFsrcIdx)),2,'omitnan');
%             RT = mean(sig(:,1:4),2,'omitnan');
%             RF = mean(sig(:,5:8),2,'omitnan'); 
%             LT = mean(sig(:,9:12),2,'omitnan'); 
%             LF = mean(sig(:,13:16),2,'omitnan'); 

            % Store the clustered ROIs
            HbORoi_t{iii}{jjj} = [RT RF LT LF];
            if jjj == 1
                pidList = [pidList; pid{iii}];
                pidAcceptIdx = [pidAcceptIdx; iii];
            end
        end

%         %For plotting waveforms for individual participants
%         subplot(2,2,jjj);
%         t = [0:1/10:15];
%         plot(t,HbORoi{iii}{jjj});
%         title(['Condition',num2str(jjj)]);
%         legend('RT','RF','LT','LF')
%         grid on;
    end
end
HbORoi = HbORoi_t(~cellfun('isempty',HbORoi_t));

%% The stats stuff: mean, standard error and plotting
for conditions = 1:4 %Loop through all four conditions
    for subjs = 1:size(HbORoi,2) %Loop through all participants
        xx(:,:,subjs) = HbORoi{subjs}{conditions}; % size: [time pts x condition x #participants]
    end

    %The average signal collapsed across all participants for each ROI for
    %every condition
    grandAvg{conditions} = mean(xx,3); 

    %The average oxygenation across all participants. Note: may want to
    %adjust size of analysis window. Do it here.
    avgOxygenation{conditions} = mean(mean(xx,1),3);

    %Calculate the standard errors
    SE{conditions} = std(mean(xx,1),0,3)/sqrt(size(HbORoi,2));
end

% Reassign variable names for readability
%   Each cond variable has 4 columns. The columns represent the following
%   regions in order: RT, RF, LT, LF
cond1 = avgOxygenation{1, 1};
SEcond1 = SE{1, 1};
cond2 = avgOxygenation{1, 2};
SEcond2 = SE{1, 2};
cond3 = avgOxygenation{1, 3};
SEcond3 = SE{1, 3};
cond4 = avgOxygenation{1, 4};
SEcond4 = SE{1, 4};

% Visualizations
fig = figure;
% Right Temporal Cluster
sph(1) = subplot(1,4,4);
RT_high = [cond1(1) cond2(1)];
RTSE_high = [SEcond1(1) SEcond2(1)];
RT_low = [cond3(1) cond4(1)];
RTSE_low = [SEcond3(1) SEcond4(1)];
errorbar(RT_high,RTSE_high,'b','LineWidth',1.5)
hold on
errorbar(RT_low,RTSE_low,'r','LineWidth',1.5)
grid on
title('Right Temporal')
xlim([0.5 2.5])
xticklabels({'Easy','Hard'})
lgd = legend('high','low','FontSize',13,'Location','southeast');
lgd.Title.String = "Context";

% Right Frontal Cluster
sph(2) = subplot(1,4,3);
RF_high = [cond1(2) cond2(2)];
RFSE_high = [SEcond1(2) SEcond2(2)];
RF_low = [cond3(2) cond4(2)];
RFSE_low = [SEcond3(2) SEcond4(2)];
errorbar(RF_high,RFSE_high,'b','LineWidth',1.5)
hold on
errorbar(RF_low,RFSE_low,'r','LineWidth',1.5)
title('Right Frontal')
xlim([0.5 2.5])
xticklabels({'Easy','Hard'})
grid on

% Left Frontal Cluster
sph(3) = subplot(1,4,1);
LF_high = [cond1(4) cond2(4)];
LFSE_high = [SEcond1(4) SEcond2(4)];
LF_low = [cond3(4) cond4(4)];
LFSE_low = [SEcond3(4) SEcond4(4)];
errorbar(LF_high,LFSE_high,'b','LineWidth',1.5)
hold on
errorbar(LF_low,LFSE_low,'r','LineWidth',1.5)
title('Left Frontal')
xlim([0.5 2.5])
xticklabels({'Easy','Hard'})
grid on

% Left Temporal Cluster
sph(4) = subplot(1,4,2);
LT_high = [cond1(3) cond2(3)];
LTSE_high = [SEcond1(3) SEcond2(3)];
LT_low = [cond3(3) cond4(3)];
LTSE_low = [SEcond3(3) SEcond4(3)];
errorbar(LT_high,LTSE_high,'b','LineWidth',1.5)
hold on
errorbar(LT_low,LTSE_low,'r','LineWidth',1.5)
title('Left Temporal')
xlim([0.5 2.5])
xticklabels({'Easy','Hard'})
grid on

linkaxes(sph,'y')
% set(sph,'YTick',[-1e-5:1e-6:5e-6])
%Give common xlabel and ylabel
han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Oxygenation (\muM.mm)','FontSize',15,'FontWeight','bold');
xlabel(han,'Signal-to-Noise Ratio','FontSize',15,'FontWeight','bold');
