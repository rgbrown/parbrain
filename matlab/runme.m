fmet = @(t, x, y) ones(size(x));
fp0 = @(t) ones(size(t));


[S, fode, u0] = treesim('N', 13);
%%
k = 2;
clear foo
Jb = sparse(ones(4*k));
for i = (1:(S.n - S.m) / k)
    foo{i} = Jb;
end
J = blkdiag(foo{:});
opts = odeset('JPattern', J);



[T, Y] = ode15s(fode, [0 500], u0, opts);
u0 = Y(end, :).';
tic
[T, Y] = ode15s(fode, linspace(0, 200, 500), u0, opts);
toc
plot(T, Y);
%%
na = S.n - S.m;
idx = S.i_spatial;
x = S.X(1:na); x = x(idx);
y = S.Y(1:na); y = y(idx);
ii = 3:4:size(Y, 2);
ngrid = sqrt(na);
Xg = reshape(x, ngrid, ngrid);
Yg = reshape(y, ngrid, ngrid);
Z = reshape(Y(1, ii(idx)), ngrid, ngrid);
h = surf(Xg, Yg, Z, 'edgecolor', 'none');
set(gca, 'zlim', [2.5 4.5])
caxis([2.5 4.5])
colormap hot
set(gcf, 'renderer', 'opengl')
%% ca
for i = 1:numel(T)
    set(h, 'ZData', reshape(Y(i, ii(idx)), ngrid, ngrid));
    drawnow
    
    
end


