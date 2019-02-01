function dic = findMultipleOutputs(filename,partialVariableNames,removeLastPoint)

load(filename)          % name of the saved file
%load T1withSST.mat

name2 = name';               % twist for right spelling
name2 = cellstr(name2(:,:)); % convert char name2 to cells

inputIndex = [];
for i=1:length(partialVariableNames)
    % find string with keyword
    a = strfind(name2,partialVariableNames{i});
    % find corresponding index
    ind = find(~cellfun(@isempty,a));
    if isempty(ind)
        str = sprintf('Variable %s not found in %s', partialVariableNames{i} , filename);
        error(str)
    end
    inputIndex = [inputIndex; ind];
end;

data = cell(length(inputIndex),1);
for i=1:length(inputIndex)
    d = findOutput(filename, name2(inputIndex(i)));
    if removeLastPoint
        data{i} = d(1:end-1);
    else
        data{i} = d;
    end;
end;
dic = containers.Map(name2(inputIndex), data);
end