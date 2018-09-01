% runPuffClassification(data, classifier, varargin) calls puffapy.py to convert MATLAB tracks struct into Python arrays,
% builds a random forest classifier off the training data and runs on the test data.
%
% This function saves the following to the Classification folder of each dataset:
%                  processedTracks.npy : all information from processedTracks.mat converted to a NumPy array
%                  RFclassifier : the Random Forest classifier generated and used
%                  RFresults.mat : the track indices classified as puffs, nonpuffs and maybes
%                  2D, 3Dfig.jpg : 2D and 3D scatter plots of results
%
% Inputs
%                  data : list of movies, using the structure returned by loadConditionData.m
%            classifier : filename of classifier in directory of training set
%                         Classifier is used if already exists or is saved as filename if it doesn't
%
% Options ('specifier', value)
%          'IsTraining' : {true}|false. Main processedTracks file is used for training
%                'File' : Name of processedTracks file (.mat or .npy) from data.source. Default: processedTracks.mat
%                         Used by default as both the training and test set
%          'SecondFile' : Name of second processedTracks file(.mat or .npy). Default: ''
%                         Must be full path to the file.
%                         If provided, used as test set by default
%              'Fields' : String of parameters to use for classification separated by spaces. Default: ''
%                         Required if training set is .mat
%                         isPuff must always be the first field
%           'Overwrite' : true|{false}. Overwrite previous tracking result.
%
% Example: runPuffClassification(data, 'RFClassifier',
%                                'IsTraining', false,
%                                'File', 'ProcessedTracks.npy',
%                                'SecondFile', 'C:\Users\Cell\Ch1\Tracking\ProcessedTracks.mat',
%                                'Fields', 'isPuff pallAdiff pfallR2 pvp')
% Notes:
% 1) Training dataset does not need to have all tracks scored.
%    puffapy.py automatically extracts and only uses the ones that have been scored.
% 2) Required Python packages: Python 3.5, H5Py, NumPy, SciPy, SkLearn, MatPlotLib
%
% Tiffany Phan, July 2016

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
ip.addParamValue('Overwrite', true, @islogical);
ip.parse(data, classifier, varargin{:});

isTraining = ip.Results.IsTraining;
relativePath = ip.Results.RelativePath;
secondFile = ip.Results.SecondFile;
file = ip.Results.File;
overwrite = ip.Results.Overwrite;

%load(file);
nd = numel(data);
for i = 1:nd
% Finds proper path to file
    [~,~,ext] = fileparts(file);
    if ext == '.npy'
        relativePath = 'Classification';
    end
    filePath = [data(i).source relativePath filesep file];

    % Finds paths to training and testing data
    if isTraining && exist(filePath, 'file')==2
        trainPath = filePath;
        if ~isempty(secondFile) && exist(secondFile, 'file')==2
            testPath = secondFile;
        elseif ~isempty(secondFile) && ~exist(secondFile, 'file')==2
            testPath = input('\n Path to second file does not exist. Please enter correct path for testing: ', 's');
        else
            testPath = '';
        end
    elseif ~isTraining && exist(filePath, 'file')==2
        testPath = filePath;
        if isempty(secondFile)
            trainPath = input('\n Please enter full path name of tracks to be used for training: ', 's');
        elseif ~isempty(secondFile) && ~exist(secondFile, 'file')==2
            trainPath = input('\n Path to second file does not exist. Please enter correct path for training: ', 's')
        else
            trainPath = secondFile;
        end
    else
        fprintf(['\n' file ' does not exist in ' data(i).source]);
        return;
    end

    %Checks if fields are entered in appropriately
    [~,~,ext] = fileparts(trainPath);
    if ext == '.mat' & isempty(ip.Results.Fields)
        fields = input('\n Please enter parameter names to train classifier with as one string separated by spaces: ' , 's');
    elseif ext == '.npy' & ~isempty(ip.Results.Fields)
        fprintf('\n Classification will be run based on fields from .npy file of training data, not the fields entered.');
        fields = '';
    else
        fields = ip.Results.Fields;
    end

    %Runs main if classification folder does not exist for test data or overwrite is true
    if ~exist(fullfile(fileparts(fileparts(testPath)),'Classification')) || overwrite
        fprintf('Running Random Forest classification (%s)...\n',getShortPath(data(i)));
        main(trainPath, testPath, fields, classifier);
        %res = [data(i).source 'Classification' filesep 'RFresults.mat'];
        %res = load(res);
        %for i = 1:numel(res.puffs)
        %    tracks(res.puffs(i)).isPuff = 1;
        %end
        %for i = 1:numel(res.nonpuffs)
        %    tracks(res.nonpuffs(i)).isPuff = 2;
        %end
        %save(filePath,'tracks','-v7.3');
    else
        fprintf('\n Classification has already been run for ', testPath);
    end
end


function main(trainPath, testPath, fields, classifier)

%Make classification folders if they don't exist
paths = {trainPath testPath};
for i = 1:numel(paths)
    if ~isempty(paths{i})
        cPath{i} = fullfile(fileparts(fileparts(paths{i})), 'Classification');
        if ~(exist(cPath{i}, 'dir')==7)
            mkdir(cPath{i});
        end
    end
end

%Properly formats inputs for command prompt
puffapy = ['"' which('puffapy.py') '"'];
classifierDir = cPath{1};
classifierPath = ['"' classifierDir filesep classifier '"'];
trainPath = ['"' trainPath '"'];
if ~isempty(fields)
    fields = ['--fields ' fields];
end
if ~isempty(testPath)
    testPath = ['--testing ' '"' testPath '"'];
end

%Call command prompt to run puffapy.py
setenv('DYLD_LIBRARY_PATH','/usr/local/lib/python3.6/site-packages/scipy/.dylibs');
systemCommand = strjoin({'/usr/local/bin/python3' puffapy classifierPath trainPath fields testPath}, ' ');
system(systemCommand);
