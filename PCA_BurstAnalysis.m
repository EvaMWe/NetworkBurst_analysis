%% PCA-ANALYSIS
% by Eva-Maria Weiss
% Version from 11.03.2020
% according to "A tutorial on Principal Components Analysis" By  Lindsay I
% Smith (2003)

% this is a PCA for one time point: data input must be a cell array:

% STEP 1: get data
% data matrix need to have the variables within the rows and the features
% withing the colums
%
% Excample
%      |mean firing rate | bursting rate | ISI |
%well_1|x11              |x12            |x13
%well_2|x21              |x22            |x23        
%well_n|x31              |x32            |x33
%----------------------------------------- 


function finData = PCA_BurstAnalysis(varargin)

%% STEP 1: get data
variableNames = {'ID','wMBR','MBD','MIBSR'};

if nargin == 1
    data = varargin{1};
elseif nargin == 2
    pathname = varargin{1};
    filename = varargin{2};
    data = load(fullfile(pathname, filename));
else
    [file,path]=uigetfile('.mat','select PCA cell data');    
    data = load(fullfile(path, file));
end

%% STEP 2:subtract the mean from each of the data dimensions

matrix_data = cell2mat(data(2:end,2:end));
matrix_data(isnan(matrix_data)) = 0;
meanVec = mean(matrix_data,1);
stdVec = std(matrix_data,0,1);
dataAdjust = (matrix_data-meanVec)./stdVec;

%% STEP 3: Calculate Covariance matrix
C = cov(dataAdjust);

%% STEP 4: Calculate Eigenvalues and Eigenvectors
[rightEigVec, eigVal] = eig(C);
eigValVec = eigVal(eigVal ~= 0); %sorted
[~, PCnb] = Maximum_median(eigValVec(:),10,'Type','number');
% get index
eigValM = repmat(eigValVec',length(PCnb),1);
indexLog = logical(sum(eigValM == PCnb));

featureVector = rightEigVec(:,indexLog);

finData = featureVector' * dataAdjust';
end
    
