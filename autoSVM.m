% (TP):Purpose of this program: builds an SVM model for puffs based on given data and desired
% parameters for comparison

% Inputs: 
% data is any Processed Tracks struct already loaded into the workspace, generated after runPuffTrackProcessing.m
% param1 and param2a re strings of the field name e.g. 'riseR2'. They will
% be plotted on the x and y-axes respectively. 

% Outputs: 
% SVMModel: saved classification SVM model generated from the input data 
% res: results for tested data, whether a track is a puff(1) vs. nonpuff(2) using the generated SVMModel

function [SVMModel, res] = autoSVM(data, param1, param2, varargin)

%if the SVMModel is to be run on new data, new data is another ProcessedTracks
%struct alredy loaded into the workspace
ip = inputParser;
ip.addParamValue('NewData', []);
ip.parse(varargin{:});

% Filter out everything we haven't classified (ZYW)
tracksFilt = data([data.isPuff] >= 1);

% Rename start and end variables in the struct
[tracksFilt.startFrame] = tracksFilt.start;
tracksFilt=rmfield(tracksFilt,'start');
[tracksFilt.endFrame] = tracksFilt.end;
tracksFilt=rmfield(tracksFilt,'end');

% Convert it to a table...because?
tracksFiltTable = struct2table(tracksFilt);

% Fit your model
tracksWanted= [[tracksFilt.(param1)]' [tracksFilt.(param2)]'];
SVMModel = fitcsvm(tracksWanted,tracksFiltTable.isPuff);

% Graphing SVMModel
sv = SVMModel.SupportVectors;
figure
gscatter(tracksWanted(:,1),tracksWanted(:,2),tracksFiltTable.isPuff);
hold on
plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
legend('Puff','Non-Puff','Support Vector')
title('SVMModel')
xlabel(param1);
ylabel(param2);
hold off

%Run SVMModel on all of current data set or new data set
newData = ip.Results.NewData; 
if isempty(newData)
    t= [data.(param1); data.(param2)];
else 
    t = [newData.(param1); newData.(param2)];
end 
t= t';

res = predict(SVMModel,t);

%Graphing Predicted Results 
figure
gscatter(t(:,1),t(:,2),res);
hold on
legend('Puff','Non-Puff')
title('Predicted Results')
xlabel(param1);
ylabel(param2);
hold off



