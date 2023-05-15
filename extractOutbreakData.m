function datai = extractOutbreakData(data, num)

idx1 = find(contains(data.Number, num2str(num)));
idx1 = idx1(1);
idx2 = find(contains(data.Number, num2str(num+1)));
idx2 = idx2(1) - 1;

datai = data(idx1:idx2, 2:3);
datai(isnat(datai.Date), :) = [];
end