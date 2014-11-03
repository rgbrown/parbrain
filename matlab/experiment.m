S = TreeSimulator('N', 13, 'RMULT', 1.2);
% Connect an nvu thing too it
fmet = @(t, x, y) ones(size(x));
[fnvu, P] = nvu('fmet', fmet);

S.nnvu = P.nVars;
S.fnvu = fnvu;
S.fp0  = @(t) ones(size(t));

%% Let's get cracking!
% Find an initial condition
u = ones(S.nnvu * S.nBlocks, 1);
f = @(t, u) S.evaluate(t, u);
S.setjpattern();
opts = odeset('JPattern', S.JPattern);
[T, U] = ode15s(f, [0 500], u, opts);
u0 = U(end, :).';

%%

fm = @(t, x, y) 1 + exp(-(x.^2 + y.^2) / ((1000e-5)^2)) * ...
    exp(-((t - 100) / 100).^2);
[fnvu, P] = nvu('fmet', fm)
S.fnvu = fnvu;
[T, U] = ode15s(f, linspace(0, 300, 200), u0, opts);
%%
if 0
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
end
%%
Pt = zeros(numel(T), S.nBlocks);
Q = zeros(numel(T), S.nBlocks);
% Compute the flow for each time step
for i = 1:numel(T)
    % Set conductance
    S.tree.setconductance(U(i, 1:S.nnvu:end).^4, 1:S.nBlocks);
    S.tree.solve(S.fp0(T(i)), S.pcap);
    Q(i, :) = S.tree.q(1:S.nBlocks);
    
end
Q = Q(:, S.tree.i_spatial(1:S.nBlocks));
disp([min(Q(:)), max(Q(:))]);

ng = sqrt(S.nBlocks);

h = imagesc(reshape(Q(1, :), ng, []), [0.14, 0.23]);
for i = 1:numel(T)
    set(h, 'Cdata', reshape(Q(i, :).', ng, []))
    pause(0.02);
    drawnow
    
end