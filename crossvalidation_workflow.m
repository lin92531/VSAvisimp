% This script performs a 12-fold cross-validation on single participant data sets 
% the statistics and machine learning toolbox and the third-party
% erpcca toolbox available at https://gitlab.com/christoph.reichert/erpcca
% are required to run the script
% data can be downloaded here: https://zenodo.org/records/TBA

experimentID = 2; % 1 = degraded vision experiment (conditionIDs 0,1,2)
                  % 2 = binocular/monocular vision experiment (conditionIDs 3,4,5)

conditionIDs = [0,1,2,4,5,6]; % list of conditionIDs to be involved
            % 0 = normal vision, 1 = MDV, 2=HDV, 3=BV, 4=MREV, 5 = MLEV

datapath = '../data/';    % path to the *.mat files

Alphabet = 'YN'; % targets were Yes or No
Contrast = true; % contrast targets and nontargets (difference wave)
ComponentLimit = 3; % limit the number of components to 3
Classifier = 'maxcorr'; % maximum correlation to determine tagets
FeatureSpace = 'correlation'; % 'maxcorr' classifier requires 'correlation' 
ReferenceType = 'impulse';

% channel list excluding TP7/TP8
chanList = {'Fp1','Fp2','F7','F3','Fz','F4','F8','FC1','FC2',...
                             'C3','Cz','C4',...
                      'PO9','CP1','CP2','PO10',...
                        'P7','P3','Pz','P4','P8',...
                      'PO7','PO3','Oz','PO4','PO8',...
                                  'IZ'};

% get list of files for all subjects
flist = dir([ datapath, 'P*.mat']);  

if isempty(flist)
    disp('Data not found. Please define the correct data path.');
end

accuracy = zeros(size( flist));

for sk = 1:length(flist)
    % load data structures
    dat = load( [ datapath, flist(sk).name]);

    % preprocess data: segment, filter and downsample data; define labels
    % and targets
    [ X,Y,T,S,meta] = preprocess_vsavisimp( dat.bciexp{experimentID},...
                                'resampfreq', 50, 'chanList', chanList);
    
    % select conditions and exclude passive vision trials
    incl = ismember(meta.cond, conditionIDs) & meta.isAttentionTrial;

    % perform 12fold cross-validation with folds comprising randomly drawn samples
    cv = crossvalerpcca( X(:,:,incl), Y(incl), T(incl), S(incl), ...
        '12-fold-rand',...
        'Alphabet', Alphabet,...
        'Contrast', Contrast,...
        'ComponentLimit',ComponentLimit,...
        'ReferenceType',ReferenceType, ...        
        'Classifier', Classifier,...
        'FeatureSpace', FeatureSpace);
    accuracy(sk) = cv.accuracy;    
end
fprintf( 'Average accuracy is %3.1f%% (std: %2.1f%%)\n',mean(accuracy*100),std(accuracy*100));
