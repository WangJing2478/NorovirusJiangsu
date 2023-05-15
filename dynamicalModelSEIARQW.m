classdef dynamicalModelSEIARQW

    properties
        N double {mustBeNonnegative, mustBeFinite}

        beta double {mustBeNonnegative, mustBeFinite}
        kappa double {mustBeNonnegative, mustBeFinite}
        p double {mustBeNonnegative, mustBeFinite}
        omega double {mustBeNonnegative, mustBeFinite}
        omegap double {mustBeNonnegative, mustBeFinite}
        gamma double {mustBeNonnegative, mustBeFinite}
        gammap double {mustBeNonnegative, mustBeFinite}
        gammapp double {mustBeNonnegative, mustBeFinite}
        q double {mustBeNonnegative, mustBeFinite}

        mu double {mustBeNonnegative, mustBeFinite}
        epsilon double {mustBeNonnegative, mustBeFinite}
        

        b
        bW

        c
        i0
        e0
        w0

        tSpan 
        xInit double {mustBeNonnegative, mustBeFinite}
        xFinal 
        t
        x

        Reff
    end

    methods
        function dxdt = derivatives(md, t, x)
            % compute the modeled derivatives for variable x at time t
            % variable vector x = [s, e, i, a, r, g, w]
            s = x(1);
            e = x(2);
            i = x(3);
            a = x(4);
            r = x(5);
            g = x(6);
            w = x(7);


            dsdt = -md.b * s * (i + md.kappa * a) - md.bW * s * w;
            dedt = -dsdt - (md.p * md.omegap + (1-md.p) * md.omega) * e;
            didt = (1-md.p) * md.omega * e - (md.q * md.gammapp + (1-md.q) * md.gamma) * i;
            dadt = md.p * md.omegap * e - md.gammap * a;
            drdt = (1-md.q) * md.gamma * i + md.gamma * g + md.gammap * a;
            dgdt = md.q * md.gammapp * i - md.gamma * g;
            dwdt = (i + a * md.c) * md.epsilon - md.epsilon * w;

            dxdt = [dsdt; dedt; didt; dadt; drdt; dgdt; dwdt];
        end

        function [t, x] = predict(md)
            % predict model with initial value md.xInit over time span md.tSpan
            odefun = @(t,x) derivatives(md, t, x);
            [t, x] = ode23(odefun, md.tSpan, md.xInit(:));
        end

        function loss = objectivefunction(parametersToBeFitted, md, tData, iData)
            % setup model with tries of parametersToBeFitted
            md.b = parametersToBeFitted(1); % parameter b
            md.bW = parametersToBeFitted(2); % parameter bW 
            md.c = parametersToBeFitted(3); % parameter c
            md.xInit(2) = parametersToBeFitted(4); % variable e
            md.xInit(3) = parametersToBeFitted(5); % variable i
            md.xInit(7) = parametersToBeFitted(6); % variable w

            % simulate (shooting)
            [t, x] = predict(md);
            dailyIncidentRate = (1-md.p) .* md.omega .* x(:,2);
            predictedCummulateI = cumtrapz(t, dailyIncidentRate);
            predictedCummulateI = interp1(t, predictedCummulateI, tData);
            observedCummulateI = cumtrapz(tData, iData);

            % losses 
            loss = norm(predictedCummulateI - observedCummulateI) ^ 2;
        end
        
        function [fittedModel, t, x] = fit(md, tData, iData, parameterNames, initialGuessOfParameters, lowerBoundOfParameters, upperBoundOfParameters)

             % fit for whole time segment
            % parameterNames = {'b', 'bW', 'c', 'e0', 'i0', 'w0'};
            if ~isempty(md.b),  b = md.b;  else b = 0;  end
            if ~isempty(md.bW),  bW = md.bW;    else bW = 0;    end
            if ~isempty(md.c),  c = md.c;  else c = 0;  end
            if ~isempty(md.e0),  e0 = md.e0;    else e0 = md.xInit(2);  end
            if ~isempty(md.i0),  i0 = md.i0;    else i0 = md.xInit(3);  end
            if ~isempty(md.w0),  w0 = md.w0;    else w0 = 0;  end

            defaultParameterValues = [b, bW, c, e0, i0, w0];
            md.tSpan = [min(tData), max(tData)];

            % idx = contains({'b', 'bW', 'c', 'e0', 'i0', 'w0'}, parameterNames);
            temp = string({'b', 'bW', 'c', 'e0', 'i0', 'w0'}) == string(parameterNames(:));
            temp = mat2cell(temp, ones(1,numel(parameterNames)), 6);
            idx = cellfun(@(x)find(x), temp);
                
            defaultParameterValues(idx) = initialGuessOfParameters;

            S.type = '()';
            S.subs = {idx};
            objfun = @(parametersToBeFitted) objectivefunction(subsasgn(defaultParameterValues, S, parametersToBeFitted), md, tData, iData);


            initialTryOfParameters = defaultParameterValues(idx);
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            lb = lowerBoundOfParameters;
            ub = upperBoundOfParameters;
            
            fittedParameters = fmincon(objfun, initialTryOfParameters, A, b, Aeq, beq, lb, ub);
            fittedParameters = subsasgn(defaultParameterValues, S, fittedParameters);

            % returns
            fittedModel = md;
            fittedModel.b = fittedParameters(1);
            fittedModel.bW = fittedParameters(2);
            fittedModel.c = fittedParameters(3);
            fittedModel.e0 = fittedParameters(4); % fitted initial e(0)
            fittedModel.i0 = fittedParameters(5); % fitted initial i(0)
            fittedModel.w0 = fittedParameters(6); % fitted initial w(0)
            fittedModel.xInit(2) = fittedModel.e0;
            fittedModel.xInit(3) = fittedModel.i0;
            fittedModel.xInit(7) = fittedModel.w0;

            [t, x] = predict(fittedModel);
            fittedModel.t = t;
            fittedModel.x = x;
            fittedModel.xFinal = x(end,:);
        end

        function [modelList, t, x] = piecewiseFit(md, tData, iData, breakPoints, parameterNames, initialGuessOfParameters, lowerBoundOfParameters, upperBoundOfParameters)

            tStart = min(tData);
            tEnd = max(tData);
            breakPoints(breakPoints <= tStart | breakPoints >= tEnd) = [];
            breakPoints = [tStart; breakPoints; tEnd];
            segmentCount = numel(breakPoints) - 1;

            modelList = repmat(md, [segmentCount, 1]);
            observedCummulation = iData(1);
            predictedCummulation = iData(1);
            t = [];
            x = [];
            for i = 1:segmentCount
                % data for this segment
                idx = tData >= breakPoints(i) & tData <= breakPoints(i+1);
                tDatai = tData(idx);
                iDatai = iData(idx);

                % init values for this segment
                model0 = md;
                model0.tSpan = [breakPoints(i), breakPoints(i+1)];
                if i > 1
                    model0.xInit = modelList(i-1).xFinal(1:7);
                    model0.e0 = modelList(i-1).xFinal(2);
                    model0.i0 = modelList(i-1).xFinal(3);
                    model0.w0 = modelList(i-1).xFinal(7);
                end


                if i > 1
                    idx = contains(parameterNames, {'e0', 'i0', 'w0'});
                    parameterNames(idx) = [];
                    initialGuessOfParameters(idx) = [];
                end


                [modeli, ti, xi] = fit(model0, tDatai, iDatai, parameterNames, initialGuessOfParameters, lowerBoundOfParameters, upperBoundOfParameters);


                modelList(i) = modeli;
                if i == 1
                    t = [t; ti];
                    x = [x; xi];
                else
                    t = [t; ti(2:end)];
                    x = [x; xi(2:end, :)];
                end
            end
        end

        function [cumulativeIncidence] = cumulativeIncidence(md)
            [t,x] = predict(md);
            cumulativeIncidence = trapz(t, (1-md.p) * md.omega * x(:, 2), 1);
        end

        function d = peak(md)
            maxDuration = 3e3;
            md.tSpan = [md.tSpan(1), maxDuration];
            [t,x] = predict(md);

            obj = (1-md.p) * md.omega * x(:, 2); % dailyNewIncidence
            d = find(obj == max(obj));
            d = t(d(end));
            if isempty(d) || d == 0 || d == maxDuration
                d = inf;
            end
        end

        function d = last(md)
            % compute the time instance of the end of epidemic (the date that daily new incidence lesser than 1)
            maxDuration = 3e3;
            md.tSpan = [md.tSpan(1), maxDuration];
            [t,x] = predict(md);

            obj = (1-md.p) * md.omega * x(:, 2); % dailyNewIncidence
            d = find(obj > 1 / md.N);

            if numel(d) == 0
                d = inf;
                return;
            end
            d = d(end);
            d = t(d);

            if isempty(d) || d == 0 || d == maxDuration 
                d = inf;
            end
        end

        function R = get.Reff(md)
        % TODO

        end
        

    end % methods

end % classdef