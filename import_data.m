%% Import data from H-tree experiment
%% Initialize variables.
filename = 'info.dat';
delimiter = ' ';
startRow = 2;
formatSpec = '%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

%% Allocate imported array to column variable names
n_processors = dataArray{:, 1};
n_blocks = dataArray{:, 2};
eqs_per_block = dataArray{:, 3};
m_local = dataArray{:, 4};
n_local = dataArray{:, 5};

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Now, import the rest
for i = 1:n_processors
    fp(i) = fopen(sprintf('out%d.dat', i-1), 'r');
end
for i = 1:n_processors
    x{i} = fread(fp(i), n_blocks, 'double');
    y{i} = fread(fp(i), n_blocks, 'double');
end    

nn = n_blocks * eqs_per_block + 1;
count = repmat(nn, n_processors, 1);
Y = repmat({zeros(nn - 1, 0)}, 1, n_processors);
t = [];
while all(count == nn)
    for i = 1:n_processors
        [s, count(i)] = fread(fp(i), nn, 'double');
        if count(i) == nn;
            if i == 1
                t = [t, s(1)];
            end
            Y{i} = [Y{i} s(2:end)];
        end
    end
    
end

%% close files
for i = 1:n_processors
    fclose(fp(i));    
end
    



