%% tabular results of piecewise-fit (write file "piecewiseFitRecord.xlsx" to current folder )
close all; clear; clc;

% 0 - not plot 1 - plot
doPlot = 1;

% disinfection (set bW = 0) in multi-routine transmission
disinfection = 0;

% initialize record
Record = cell(0);


%% Read data
opts1 = detectImportOptions("allData.xlsx", 'Sheet', '流行曲线');
opts2 = detectImportOptions("allData.xlsx", 'Sheet', '疫情基本信息');
data = readtable('allData.xlsx', opts1, 'Sheet', '流行曲线');
dataInfo = readtable('allData.xlsx', opts2, 'Sheet', '疫情基本信息');



%% main loop, analyze each outbreak
for outbreakNumber = 1:23


    %% initialize model
    model = dynamicalModelSEIARQW;

    % select and extract data for the outbreak to be analyzed
    model.N = dataInfo.PopulationSize(dataInfo.ID == outbreakNumber);
    datai = extractOutbreakData(data, outbreakNumber);
    datai(datai.Incidence == 0, :) = [];

    %datai = readtable('data.xlsx');
    iData = datai.Incidence / model.N;
    dData = datai.Date;
    tData = days(dData - min(dData));

    if numel(tData) == 0
        failedFlag = "missingObservations";
        temp = cell(1,10);
        temp{1} = outbreakNumber;
        temp{2} = dataInfo.Routine(dataInfo.ID == outbreakNumber);
        temp{4} = failedFlag;
        Record = [Record; temp];
        continue;
    end

    % data interpolation
    [tRefined, iRefined] = deal(tData, iData);

    %% Setup model
    % basic parameters
    variableNames = {'s', 'e', 'i', 'a', 'r', 'q', 'w'};
    model.kappa = 0.05;
    model.omega = 1;
    model.omegap = 1;
    model.p = 0.3;
    model.gamma = 1 / 3;
    model.gammap = 0.03846;
    model.gammapp = 2;
    model.epsilon = 0.1;
    model.q = 0;

    % default intial value of variables
    model.e0 = iData(1) / ((1-model.p) * model.omega);
    model.i0 = iData(1);
    model.w0 = 0;
    model.xInit = [1, model.e0, model.i0, 0, 0, model.w0];

    %% Setup optimization problem
    if strcmp(dataInfo.Routine(dataInfo.ID == outbreakNumber), '水/食物传人')
        parametersToBeFitted = {'b', 'bW', 'w0', 'c'};
        initialGuessOfParameters = [1e-1, 1e-1, 1e-1, 0.3];
        lowerBoundOfParameters = [0, 0, 1e-3, 1e-3];
        upperBoundOfParameters = [1, 1, 1, 1];
    else
        %     parametersToBeFitted = {'b','i0'};
        %     initialGuessOfParameters = [0.5, 1e-4];
        %     lowerBoundOfParameters = [0, 0];
        %     upperBoundOfParameters = [1000, 1e-4];

        parametersToBeFitted = {'b'};
        initialGuessOfParameters = [3];
        lowerBoundOfParameters = [0];
        upperBoundOfParameters = [100];

    end

    %% piecewise fit


    if numel(tData) >= 2  % at least 2 observations
        breakPoints = tData(iData == max(iData));
        breakPoints = breakPoints(1);
        [modelList, t, x] = piecewiseFit(model, tRefined, iRefined, breakPoints, parametersToBeFitted, initialGuessOfParameters, lowerBoundOfParameters, upperBoundOfParameters);
        dailyNewIncidence = (1-model.p) * model.omega * x(:,2); % piecewise fitted incidence

        if modelList(1).b == initialGuessOfParameters(1)
            failedFlag = "stationaryInit";
        else
            failedFlag = "success";
        end

    else
        failedFlag = "limitedObservations";

    end



    % save Record
    for i = 1:numel(modelList)
        temp = cell(1,10);
        temp{1} = outbreakNumber;
        temp{2} = dataInfo.Routine(dataInfo.ID == outbreakNumber);
        temp{3} = i;
        temp{4} = failedFlag;
        temp{5} = modelList(i).b;
        temp{6} = modelList(i).bW;
        temp{7} = modelList(i).w0;
        temp{8} = modelList(i).c;
        temp{9} = modelList(i).e0;
        temp{10} = modelList(i).i0;
        Record = [Record; temp];
    end



end

RecordTable = cell2table(Record, 'VariableNames', {'OutbreakNumber', 'Routine', 'SegmentNumber', 'FailedFlag', 'b', 'bW', 'w0', 'c', 'e0', 'i0'});
writetable(RecordTable, 'piecewiseFitRecord.xlsx');