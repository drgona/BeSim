function run_validation(reload,nStep,saveFig,disturbanceType)
    if nargin < 1
        clear all; clc; close all;
        reload = 1;
        nStep = 360*24;
        saveFig=1;
        disturbanceType='';
    end;
addpath('../../utilities/');
validation(reload, nStep, saveFig, disturbanceType);
end
