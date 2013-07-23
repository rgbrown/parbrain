function [S, f_next, f_close] = data_reader()

%% Get the configuration from info.dat
[n_processors,n_blocks,eqs_per_block,m_local,n_local,m_global,n_global] = ...
    loadinfo('info.dat');
%%
% Open the data files for reading
for i = 1:n_processors
    fp(i) = fopen(sprintf('out%d.dat', i-1), 'r');
end

% Read in the x and y coordinates of each block
for i = 1:n_processors
    x{i} = fread(fp(i), n_blocks, 'double');
    y{i} = fread(fp(i), n_blocks, 'double');
end
% Convert x and y to the right shape
reshape_fcn = @(x) cell2mat(...
    reshape(cellfun(@(s) reshape(s, m_local, n_local), x, ...
    'uniformoutput', false), m_global, n_global));
S.x = reshape_fcn(x);
S.y = reshape_fcn(y);
f_next = @get_next;
f_close = @close_files;
nread = n_blocks * eqs_per_block + 1;
%%
    function [t, Y] = get_next()
        % Read a line of raw data
        for i = 1:n_processors
            [s, count] = fread(fp(i), nread, 'double');
            if count ~= nread
                t = [];
                Y = [];
                return
            end
            if i == 1
                t = s(1);
            end
            for j = 1:eqs_per_block
                Y{j}{i} = s((j+1):eqs_per_block:end);
            end
            
        end
        for j = 1:eqs_per_block
            Y{j} = reshape_fcn(Y{j});
        end
        
    end


    function close_files()
        
        %% close files
        for i = 1:n_processors
            fclose(fp(i));
        end
        
        
    end
end

