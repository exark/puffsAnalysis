% runPuffClassification(data, classifier, varargin) calls puffapy.py to train and test Random Forest Classifier on datasets
% This function generates the following in the Classification folder of each dataset: 
%                  processedTracks.npy : all information from processedTracks.mat converted to a NumPy array
%                  RFclassifier : the Random Forest classifier generated (from training set) and used
%                  RFresults.mat : the track indices classified as puffs, nonpuffs and maybes
%                  2D,3Dfig.jpg : 2D and 3D (if possible) scatter plots of results               
%
% Inputs   
%                  data : list of movies, using the structure returned by loadConditionData.m
%            classifier : filename of classifier (to use if already exists, to save as if it doesn't)   
%
% Options ('specifier', value)
%          'IsTraining' : {true}|false. Main processedTracks file is used for training
%                'File' : Name of main processedTracks file (.mat or .npy). Default: processedTracks.mat from data.source
%                         Used by default as both the training and test set
%          'secondFile' : Name of second processedTracks file(.mat or .npy). Default: ''
%                         Must be full path to the file. 
%                         If provided, used as test set by default
%              'Fields' : String of parameters to use for classification separated by spaces. Default: ''
%                         Required if training set is .mat 
%           'Overwrite' : true|{false}. Overwrite previous tracking result.
%
% Example: runPuffClassification(data, 'RFClassifier', 
%                                'IsTraining', false, 
%                                'SecondFile', 'C:\Users\Cell\Ch1\Tracking\ProcessedTracks.mat', 
%                                'Fields', 'isPuff pallAdiff pfallR2 pvp')

% Francois Aguet, May 2010 (last modified 05/28/2013)

function runPuffClassification(data, classifier, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('classifier', @ischar);
ip.addParamValue('IsTraining', true, @islogical);
ip.addParamValue('File', 'ProcessedTracks.mat', @ischar); 
ip.addParamValue('SecondFile', '', @ischar);
ip.addParamValue('Fields', '', @ischar);
ip.addParamValue('RelativePath', 'Tracking', @ischar);
%ip.addParamValue('Overwrite', true, @islogical); % ADD THIS FUNCTIONALITY 
%ADD IN PARFOR

ip.parse(data, classifier, varargin{:});
isTraining = ip.Results.IsTraining; 
rPath = ip.Results.RelativePath;
secondFile = ip.Results.SecondFile; 
file = ip.Results.File;

% Set the file paths of training and test sets
[~,~,ext] = fileparts(file);
if ext == '.npy'
    rPath = 'Classification';
end 

tPath = [data.source rPath filesep file];
if isTraining && exist(tPath, 'file')==2
    trainPath = tPath;
    if ~isempty(secondFile) && exist(secondFile, 'file')==2
        testpath = ['--testing ' '"' secondFile '"'];
    else
        testPath = secondFile;
    end 
elseif ~isTraining && exist(tPath, 'file')==2
    testPath = ['--testing ' '"' tPath '"'];
    if isempty(secondFile) 
        trainPath = input('Please enter full path name of tracks to be used for training: ', 's');
    elseif ~isempty(secondFile) && ~exist(secondFile, 'file')==2
        fprintf([secondFile ' does not exist.']);
        trainPath = input('Please enter correct full path name of tracks to be used for training: ', 's')
    else
        trainPath = secondFile;
    end 
else 
    fprintf([file ' does not exist in that directory']);
    return; 
end

%Checks if fields are entered in appropriately 
[~,~,ext] = fileparts(trainPath);
if ext == '.mat' & isempty(ip.Results.Fields)
    fields = input('\n Please enter parameter names to train classifier with, separated by spaces: ' , 's');
    fields = ['--fields ' fields];
elseif ext == '.npy' & ~isempty(ip.Results.Fields)
    fprintf('\n Classification automatically run based on fields in the training .npy file, not the fields entered.');
    fields = ''; 
else
    fields = ip.Results.Fields; 
    fields = ['--fields ' fields];
end 

%Make classification folder if it doesn't already exist
if ~(exist([data.source 'Classification'], 'dir')==7)
    mkdir([data.source 'Classification']);
end

classifierDir = data.source;
if ~isTraining 
    classifierDir = fileparts(fileparts(secondFile));
    if ~(exist([classifierDir filesep 'Classification'], 'dir')==7)
        mkdir([classifierDir filesep 'Classification']);
    end
end

%Properly formats inputs for command prompt
classifierPath = ['"' classifierDir filesep 'Classification' filesep classifier '"'];
trainPath = ['"' trainPath '"'];
puffapy = ['"' which('puffapy.py') '"'];

%Call command prompt to run puffapy.py
systemCommand = strjoin({'python' puffapy classifierPath trainPath fields testPath}, ' ');
fprintf('\n Random Forest Classification is running...');
system(systemCommand);
