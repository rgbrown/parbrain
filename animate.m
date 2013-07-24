[S, fn, fc] = data_reader();
[t, Y] = fn();
idx = 4;
%%
h = imagesc(Y{idx}, [0 10]);


%%
while 1
    tprev = t;
    [t, Y] = fn();
    if isempty(Y)
        break
    end
    set(h, 'CData', Y{idx});
    pause(t - tprev);
end
fc();