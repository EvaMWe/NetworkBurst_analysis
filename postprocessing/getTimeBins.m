function timepoints = getTimeBins(varargin)

if nargin ~= 0
    timepoints = varargin {1};
else
    prompt = {sprintf('Enter starts for binning')};
    starts_ = inputdlg(prompt);
    starts = str2double(strsplit(starts_{1,1},','));
    
    prompt = {sprintf('Enter  stops for binning')};
    stops_ = inputdlg(prompt);
    stops = str2double(strsplit(stops_{1,1},','));
    
    timepoints =[starts; stops]';
end
end
