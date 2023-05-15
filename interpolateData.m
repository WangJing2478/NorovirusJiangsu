function [xRefined, yRefined] = interpolateData(x, y, stepSize, method)
xRefined = min(x) : stepSize : max(x);
yRefined = interp1(x, cumsum(y), xRefined, method);
yRefined = [y(1), diff(yRefined) / stepSize];
end