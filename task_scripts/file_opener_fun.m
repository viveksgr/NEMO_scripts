function chunk_x = file_opener_fun(chunk,path_,delimiter)
% chunk is a whole number: 0,1,2,3...
if nargin<3
    delimiter = ',';
end

if nargin<2
    path_ = pwd;
end
filename = [path_ '\chunk'...
            num2str(chunk) '.csv'];
startRow = 2;
formatSpec = '%f%s%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
dataArray([1, 3, 4]) = cellfun(@(x) num2cell(x), dataArray([1, 3, 4]), 'UniformOutput', false);
chunk_x = [dataArray{1:end-1}];
end

% eval(sprintf('', chunk));
% chunk_x = eval(sprintf('chunk%d', chunk));

