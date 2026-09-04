function [X,Y,T,S,meta] = preprocess_vsavisimp(bciexp, varargin)

% define parameters
prm = struct;
dsf = 1; % downsampling factor
prm.passband = [1, 12.5]; % cutoff frequencies for bandpass filter 
prm.stopband = [ 48 52];
prm.interval_length = 0.75; % length of analysis interval in seconds
prm.baseline_start = 0;
prm.baseline_length = 0;
prm.interval_start = 0;

for k=1:2:length(varargin)
    prm.(varargin{k}) = varargin{k+1};
end

if isfield(prm,'resampfreq')
    dsf = bciexp.srate/prm.resampfreq;
else
    prm.resampfreq = bciexp.srate/dsf;
end
if isfield(prm,'chanList') && ~isempty(prm.chanList)
    chanIdx = zeros(size(prm.chanList));
    for k=1:length(chanIdx)
        pos = find( strcmp(prm.chanList{k}, bciexp.label));
        if isempty(pos)
            error('Channel name not found.');
        end
        chanIdx(k)=pos;
    end
    chanIdx = sort(chanIdx);
else
    chanIdx = [1:28,30]; 
end

% rereference
lmast = find(strcmp(bciexp.label,'LMAST'));
if isempty(lmast)
    warning('No LMAST found for rereferencing.')
    return;
end
curIdx = setdiff(1:size(bciexp.data,1),lmast);
bciexp.data(curIdx,:,:) = bciexp.data( curIdx,:,:) - bciexp.data(lmast,:,:)/2;

% bandpass filter coefficients
[butterB, butterA ] = butter( 2, prm.passband / ( bciexp.srate / 2)); 

% notch filter coefficients
[notchB,notchA] = butter( 2, prm.stopband / ( bciexp.srate/2),'stop');

% determine number of samples for baseline window
nBlSamp = floor( prm.baseline_length * bciexp.srate/dsf)-1;
startSampBl = round(prm.baseline_start * bciexp.srate/dsf);

% define baseline interval 
baselineInterval = startSampBl + (0:nBlSamp-1);

% determine number of samples per analysis window
nWinSamp = floor( prm.interval_length * bciexp.srate/dsf)-1;
startSamp = round(prm.interval_start * bciexp.srate/dsf);

% define interval relative to stimulus onset
interval = startSamp + (0:nWinSamp-1);

% number of stimuli per trial
nStimPerTrial = sum( diff( sum( bciexp.stim( :, :, 1))) > 0);

% number of available trials
nTrials = size( bciexp.data,3); 

% number of channels
nChan = length( chanIdx); 

% Initialize the matrices for  segmented data and model functions
nStim = nStimPerTrial*nTrials;
X = zeros(nChan,nWinSamp,nStim);
Y = zeros(nStim,1);
T = zeros(nStim,1);
S = zeros( nStim,1);

meta.time = linspace(1/prm.resampfreq,prm.interval_length-1/prm.resampfreq,size(X,2));
meta.trialOnsets = zeros( nTrials,1);

meta.cond = zeros(nStim,1);
meta.runID = zeros(nStim,1);
meta.sessID = zeros(nStim,1);

% loop over all trials
for tr=1:nTrials
    
    % filter the data

    fdat = filtfilt( notchB, notchA, bciexp.data( chanIdx, :, tr)');

    % band pass filter
    fdat = filtfilt( butterB, butterA, fdat);
    %resample
    if dsf>1
        rdat = resample( fdat, 1, dsf);
    else
        rdat = fdat;
    end
        
    % determine stimulus onsets
    stimOnsets = find( sum(bciexp.stim(:,:,tr))~=0 );
    meta.trialOnsets(tr) = (tr-1)* length( stimOnsets)+1;
    if length( stimOnsets)~=nStimPerTrial
        error( 'Something is wrong with number of stim onsets.');
    end
    % loop over all stimuli in a trial
    for st = 1: length( stimOnsets)        
        cnt = st+(tr-1)*nStimPerTrial;
        % determine interval
        idx = ceil( stimOnsets( st) / dsf) + interval;
        
        % cut out the resampled data
        X(:,:,cnt) = rdat(idx,:)';

        % do baseline correction
        if ~isempty(baselineInterval)
            bidx = ceil( stimOnsets( st) / dsf) + baselineInterval;
            X(:,:,cnt) = X(:,:,cnt) - mean(rdat(bidx,:),1)';
        end
        % store the conditions        
        Y(cnt) = 2-bciexp.stim(1,stimOnsets(st),tr); % 1= green left 2=red left
        switch bciexp.expected{tr}
            case 'yes',  T(cnt) = 1; 
            case 'no', T(cnt) = 2;
            otherwise % '' = idle, T(cnt) remains 0 
        end
        S(cnt) = tr;
        meta.cond(cnt) = bciexp.cond(tr);
        meta.runID(cnt) = bciexp.runID(tr);
        meta.sessID(cnt) = bciexp.sessID(tr);
    end
end

meta.isAttentionTrial=T>0;