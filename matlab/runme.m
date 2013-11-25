fmet = @(t, x, y) ones(size(x));
fp0 = @(t) ones(size(t));


[S, fode, u0] = treesim('N', 9);
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
[T, Y] = ode15s(fode, [0 500], u0, opts);
plot(T, Y);

na = S.n - S.m;
idx = S.i_spatial;
x = S.X(1:na); x = x(idx);
y = S.Y(1:na); y = y(idx);
