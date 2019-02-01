function [t, TLow, TUp, TSup, TRefControl] = comfortTemperature(pathData)
    % Comfort temperature based on Norm ISO7730 and supply temperature form
    % heat curve.
%     pathData = 'Data/';
    
    fileName = [pathData 'preComputed.mat'];
    varNames = {'TZoneSetLow','TZoneSetUp','TSup','TRefControl'};

    t = findOutput(fileName, 'Time' )';
    
    try
        TLow = findOutput(fileName,varNames(1))';
        TUp = findOutput(fileName,varNames(2))';
        TSup = findOutput(fileName,varNames(3))';
        TRefControl = findOutput(fileName,varNames(4))';
    catch
        warning('Variables for comfort zone not found in ../preComputed.mat');
        TLow = zeros(size(t));
        TUp = zeros(size(t));
        TSup = zeros(size(t));
        TRefControl = zeros(size(t));
    end
end