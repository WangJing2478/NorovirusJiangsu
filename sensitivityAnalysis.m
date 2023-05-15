%% Sensitivity Analysis
clear; close all; clc;

% disinfection (set bW = 0) in multi-routine transmission
disinfection = 0;

% select one outbreak to analyze
outbreakNumber = 1;

% select parameter to be analyzed
%targetParameter = 'b';
%targetedValues = linspace(1e-3,1e0, 10);
targetParameter = 'epsilon';
targetedValues = linspace(0.03571,0.1429, 10);

%% Read data
opts1 = detectImportOptions("allData.xlsx", 'Sheet', '流行曲线');
opts2 = detectImportOptions("allData.xlsx", 'Sheet', '疫情基本信息');
data = readtable('allData.xlsx', opts1, 'Sheet', '流行曲线');
dataInfo = readtable('allData.xlsx', opts2, 'Sheet', '疫情基本信息');


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



%% setup breakPoints
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




%% initialize figure
fig3 = figure;
fig3.WindowState = 'maximized';
tiled3 = tiledlayout('flow');

%% fit for first segment (assumed without any intervention)
[model1, t1, x1] = fit(model, tRefined(tRefined <= breakPoints(1)), iRefined(tRefined <= breakPoints(1)), parametersToBeFitted, initialGuessOfParameters, lowerBoundOfParameters, upperBoundOfParameters);
dailyNewIncidence1 = (1-model.p) *  model1.omega * x1(:,2);

nexttile;
hold on;
plot(tData, iData, 'bo','LineWidth',2);
plot(t1, dailyNewIncidence1,'r-','LineWidth',2);

lastRecord = zeros(numel(targetedValues)+1,1)
for i = 1:numel(targetedValues)+1
    % i==1 for non-intervention

    model2 = model1;
    model2.xInit = model1.x(end, 1:7);
%     model2.xInit = model1.x(1, 1:7);
    %model2.tSpan = [breakPoints(1), max(tRefined)];
   
    S.type = '.';
    S.subs = targetParameter;
    if i > 1
        %model2.bW = targetedValues(i-1);
        model2 = subsasgn(model2, S, targetedValues(i-1));
    if strcmp(targetParameter,'omega')
        model2.xInit(2) = model2.xInit(2)*model1.omega/targetedValues(i-1);    
    end
    end

    model2.tSpan = [breakPoints(1), 80];
%     model2.tSpan = [0, last(model1)];
    %model2.tSpan = [0, last(model2)];
  lastRecord(i)=model2.tSpan(2);

    [t2, x2] = predict(model2);
    dailyNewIncidence2 = (1-model.p) * model2.omega * x2(:,2);

    %% series of plots
    if i == 1
        plot(t2, dailyNewIncidence2, 'black-.','LineWidth',2);
    else
        plot(t2, dailyNewIncidence2, 'r-.','LineWidth',2);
    end
ylim([0,0.05])
end

%scatter the end of epidemics
%plot(lastRecord,zeros(numel(lastRecord,1)),'blackp')

xlabel(gca, 'Time (in days)');
ylabel(gca, 'Daily New Incidence');
%legend(["Observations", "Fitted Curve", model1.bW, targetedValues], 'Location', 'northwest');
legend(["Observations", "Fitted Curve", subsref(model1, S), targetedValues], 'Location', 'northwest');
%title(tiled3, "Sensitivity Analysis for parameter " + string(targetParameter) + " (outbreak " + outbreakNumber + ", " + dataInfo.Routine(dataInfo.ID == outbreakNumber) + ")" );
title(tiled3, "Sensitivity Analysis for parameter " + string(targetParameter))

temp = [[t1; t2], [dailyNewIncidence1; dailyNewIncidence2]];
            tempTable = array2table(temp, 'VariableNames', {'Time', 'DailyNewIncidence'});
            writematrix(temp, 'sensitivity.xlsx');
            writetable(tempTable, 'sensitivity.xlsx');