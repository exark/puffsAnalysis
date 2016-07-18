% runPuffClassification(data, classifier, varargin) calls puffapy.py to train and test Random Forest Classifier on datasets
% This function generates the following in the Classification folder of each dataset: 
%                  processedTracks.npy : all information from processedTracks.mat converted to a NumPy array
%                  RFclassifier : the Random Forest classifier generated and used
%                  RFresults.mat : the track indices (corresponds to processedTracks.mat) classified as puffs, nonpuffs and maybes
%                  2D,3Dfig.jpg : 2D and 3D (if possible) scatter plots of results               
%
% Inputs   
%                  data : list of movies, using the structure returned by loadConditionData.m
%            classifier : filename of classifier (to use if already exists, to save as if it doesn't)   
%
% Options ('specifier', value)
%          'IsTraining' : {true}|false. Tracks file being passed in is used for training
%           'File'      : Name of tracks file present in currentdir(.mat or .npy). Default: ProcessedTracks.mat
%           'OtherFile' : Name of other tracks file(.mat or .npy). Default: ''
%                         Used as test set if isTraining == false. 
%                         Must be full path to the file. 
%              'Fields' : String of parameters to use for classification separated by spaces. Default: ''
%                         Required if training set is .mat 
%           'Overwrite' : true|{false}. Overwrite previous tracking result.
%
% Example: runPuffClassification(data, 'RFClassifier', 'IsTraining', false, 'TrainData', 'C:\User\test.mat',  ) ;

% Francois Aguet, May 2010 (last modified 05/28/2013)

function runPuffClassification(data, classifier, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('classifier', @ischar);
ip.addParamValue('IsTraining', true, @islogical);
ip.addParamValue('File', 'ProcessedTracks.mat', @ischar); 
ip.addParamValue('OtherFile', '', @ischar);
ip.addParamValue('Fields', '', @ischar);
ip.addParamValue('RelativePath', 'Tracking', @ischar);
ip.addParamValue('Overwrite', true, @islogical); % ADD THIS FUNCTIONALITY 
%ADD IN PARFOR

ip.parse(data, classifier, varargin{:});
isTraining = ip.Results.IsTraining; 
rPath = ip.Results.RelativePath;
otherFile = ip.Results.OtherFile; 
file = ip.Results.File;

% Set the file paths of training and test sets
[~,~,ext] = fileparts(file);
if ext == '.npy'
    rPath = 'Classification';
end 

tPath = [data.source filesep rPath filesep file];
if isTraining && exist(tPath, 'file')==2
    trainPath = tPath;
    if ~isempty(otherFile) && exist(otherFile, 'file')==2
        testpath = ['--testing ' '"' otherFile '"'];
    else
        testPath = otherFile;
    end 
elseif ~isTraining && exist(tPath, 'file')==2
    testPath = ['--testing ' '"' tPath '"'];
    if isempty(otherFile) 
        trainPath = input('Please enter full path name of tracks to be used for training: ', 's');
    elseif ~isempty(otherFile) && ~exist(otherFile, 'file')==2
        fprintf([otherFile ' does not exist.']);
        trainPath = input('Please enter correct full path name of tracks to be used for training: ', 's')
    else
        trainPath = otherFile;
    end 
else 
    fprintf([file ' does not exist in that directory']);
    return; 
end

%Checks if fields are entered in appropriately 
[~,~,ext] = fileparts(trainPath);
if ext == '.mat' & isempty(ip.Results.Fields)
    fields = input('Please enter parameter names to train classifier with, separated by spaces: ' , 's');
elseif ext == '.npy' & ~isempty(ip.Results.Fields)
    fprintf('Classification automatically run based on parameters in the training .npy file');
    fields = ''; 
else
    fields = ip.Results.Fields; 
end 

%Make classification folder if it doesn't already exist
if ~(exist([data.source 'Classification'], 'dir')==7)
    mkdir([data.source 'Classification']);
end

%Properly formats inputs for command prompt
classifierPath = ['"' data.source filesep 'Classification' filesep classifier '"'];
trainPath = ['"' trainPath '"'];
fields = ['--fields ' fields];
puffapy = ['"' which('puffapy.py') '"'];

%Call command prompt to run puffapy.py
systemCommand = strjoin({'python' puffapy classifierPath trainPath fields testPath}, ' ');
fprintf('Random Forest Classification is running...');
system(systemCommand);
