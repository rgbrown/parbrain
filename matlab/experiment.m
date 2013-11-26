S = TreeSimulator('N', 15);
% Connect an nvu thing too it
[fnvu, P] = nvu();

S.nnvu = P.nVars;
S.fnvu = fnvu;
S.fp0  = @(t) ones(size(t));

%% Let's get cracking!
% Find an initial condition
u = ones(S.nnvu * S.nBlocks, 1);
f = @(t, u) S.evaluate(t, u);
S.set_jpattern();
opts = odeset('JPattern', S.JPattern);
[T, U] = ode15s(f, [0 500], u, opts);
u0 = U(end, :).';

%%
[T, U] = ode15s(f, linspace(0, 250, 200), u0, opts);
ii = S.imask() + P.iSMC;
z = U(1, :);
z = z(ii);
figure(1)
clf
h = imagesc(z, [0 1]);
colormap hot
for i = 1:numel(T)
    z = U(i, :);
    z = z(ii);
    set(h, 'Cdata', z)
    pause(0.05);
end