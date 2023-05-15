%% fit and predict with interventions
close all; clear; clc;

% disinfection (set bW = 0) in multi-routine transmission
disinfection = 0;

% select one outbreak to analyze
outbreakNumber = 1;




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
fig2 = figure;
fig2.WindowState = 'maximized';
tiled2 = tiledlayout(1, 3);

% fit for first segment (assumed without any intervention)
[model1, t1, x1] = fit(model, tRefined(tRefined <= breakPoints(1)), iRefined(tRefined <= breakPoints(1)), parametersToBeFitted, initialGuessOfParameters, lowerBoundOfParameters, upperBoundOfParameters);
dailyNewIncidence1 = (1-model.p) *  model1.omega * x1(:,2);

predictedCumulativeIncidence = cumulativeIncidence(model1)
predictedPeak = peak(model1)
predictedLast = last(model1)

% simulate for 3 kinds of interventions
interventionNames = {'Quarantine', 'School Closure', 'Sterilization'};
for k = 1:numel(interventionNames)

    % plot data and fitted curve in the first segment
    ax = nexttile;
    hold on;
    plot(tData, iData, 'bo');
    plot(t1, dailyNewIncidence1, 'r-');

    % series of interventions in the second segment
    intensity = 0 : 0.1 : 1;
     lastRecord = zeros(numel(intensity)+1,1)
    for i = 1:numel(intensity)
        model2 = model1;
        model2.xInit = model1.x(end, 1:7);
        %model2.xInit = model1.x(1, 1:7);
        %model2.tSpan = [breakPoints(1), max(tRefined)];
        
        if model2.b >= 1
            % b >= 1, the exponential contact-probability decomposition has no real solution
            e = model2.b / 15;
            switch k
                case 1
                    % quarantine
                    model2.q = intensity(i);
                case 2
                    % school closure
                    h = 1-intensity(i);
                    model2.b = h * e * 15;
                case 3
                    % disinfecting
                    x = 1 - intensity(i);
                    model2.b = x * e * 15;
            end
        else
            % 0 < b < 1, the exponential contact-probability decomposition has unique real solution
            e = 1 - (1-model2.b) ^ (1/15);
            switch k
                case 1
                    % quarantine
                    model2.q = intensity(i);
                case 2
                    % school closure
                    h = 1-intensity(i);
                    model2.b = 1 - (1-e) ^ (15*h);
                case 3
                    % disinfecting
                    x = 1 - intensity(i);
                    model2.b = 1 - (1 - x * e) ^ 15;
            end
        end

        if disinfection == 1
                 model2.bW = 0;
        end

%model2.tSpan = [0, last(model2)];
model2.tSpan = [breakPoints(1), last(model2)];
lastRecord(i)=model2.tSpan(2);
        [t2, x2] = predict(model2);


        dailyNewIncidence2 = (1-model.p) * model2.omega * x2(:,2);

        if k == 1 && i == 11
            temp = [[t1; t2], [dailyNewIncidence1; dailyNewIncidence2]];
            tempTable = array2table(temp, 'VariableNames', {'Time', 'DailyNewIncidence'});
            writematrix(temp, 'interventionMatrix_k1i2.xlsx');
            writetable(tempTable, 'interventionTable_k1i2.xlsx');
        end 

        % series of plots
        if i == 1
            plot(t2, dailyNewIncidence2, 'black-.');
        else
            plot(t2, dailyNewIncidence2, 'r-.');
        end

    end
    %scatter the end of epidemics
    plot(lastRecord,zeros(numel(lastRecord,1)),'blackp');

    xlabel(ax, 'Time (in days)');
    ylabel(ax, 'Daily New Incidence');
    legend1 = {'Observations', 'Fitted Incidence', 'Predication Without Intervention'};
    legend2 = "Predication With Intervention Intensity " + intensity(2:end) * 100 + "%";
    legend([legend1, legend2], 'Location', 'northeast');
    title(interventionNames{k});
    ax.YLim(1) = 0;

end

title(tiled2, ["outbreak " + outbreakNumber + ", " + dataInfo.Routine(dataInfo.ID == outbreakNumber)]);
