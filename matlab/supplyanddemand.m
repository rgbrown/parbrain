%% Autoregulation function design
% We will first consider the range of tissue concentrations that can be
% produced by our model

% First create a tree - maybe 13 levels is adequate
% 1/0.6 is good for RMULT
alpha = 0.6;
S = TreeSimulator('N', 13, 'RMULT', 1/alpha, 'NS', 8);
[fnvu, params] = nvu('REGDISABLE', true, 'GAMMA', 50);
% Connect the pieces up
S.fnvu = fnvu;
S.nnvu = params.nVars;
S.fp0 = @(t) ones(size(t));

% Compute a baseline initial condition
u0 = repmat(0.5, S.nnvu * S.nBlocks, 1);
[u0, success] = S.computeeq(u0);
S.pressureflow(u0, 1);
qbase = S.tree.q(1);
cstar = u0(params.iT);
gamma = 100;
S.fnvu = nvu('REGDISABLE', true, 'CSTAR', cstar, 'GAMMA', gamma);
u0 = S.computeeq(u0);
iq = reshape(S.tree.i_spatial, sqrt(S.nBlocks), []);
X = S.x(iq);
Y = S.y(iq);
x = X(1, :);
y = Y(:, 1).';

%% Pressure regulation
% We perform an experiment where the increased for a period of time,
% allowed to return to equilibrium, and then decreased for a period of
% time. First, we simulate with regulation turned off to show what the flow
% and concentration would do
regdisabled = false;
S.fnvu = nvu('REGDISABLE', regdisabled, 'CSTAR', cstar, 'GAMMA', gamma);
a = 2;
fp0 = @(t) 1 + ...
    0.5*(1./(1 + exp(-a*(t - 20))) - 1./(1 + exp(-a*(t - 80)))) + ...
    -0.25*(1./(1 + exp(-a*(t - 120))) - 1./(1 + exp(-a*(t - 180))));
t = linspace(0, 200, 1000);
plot(t, fp0(t));
S.fp0 = fp0;

opts = odeset('JPattern', S.JPattern);
fn = @(t, u) S.evaluate(t, u);
[T, U] = ode15s(fn, [0 200], u0, opts);

%compute the flow
Q = zeros(S.nBlocks, size(U, 1));
for i = 1:size(U, 1);
   S.pressureflow(U(i, :).', fp0(T(i)));
   Q(:, i) = S.tree.q(1:S.nBlocks);
end
%%
subplot(4, 1, 1)
plot(T, fp0(T), 'k', 'linewidth', 2);
set(gca, 'fontsize', 16)
ylabel('p_0(t)')

subplot(4, 1, 2)
plot(T, Q(1, :), 'k', 'linewidth', 2);
set(gca, 'fontsize', 16)
ylabel('flow')

subplot(4, 1, 3)
plot(T, U(:, params.iT), 'k', 'linewidth', 2);
set(gca, 'fontsize', 16)
ylabel('Tissue [CO_2]')

subplot(4, 1, 4)
plot(T, U(:, params.iSMC), 'k', T, U(:, params.iR), 'k--', 'linewidth', 2);
set(gca, 'fontsize', 16, 'ylim', [0 1.5])
xlabel('t (sec)')

%%
if regdisabled
    print -deps regdisabled.eps
else
    print -deps regenabled.eps
end
%
%

end

%% Now, we'll look at changes with neurological activity levels
% Create a metabolic function that starts in one corner, and we're just
% going to plot results for a cross section along one edge

Radii = [5e-4, 1e-3, 2e-3, 3e-3, 5e-3];
nr = numel(Radii);
clear ua p q
xmax = max(x);
ymax = max(y);
for i = 1:nr

radius = Radii(i);
xmin = x(1); ymin = y(1);
fmet = @(t, x, y) 1 + 0.5*exp(-((x - xmax/2).^2 + (y - ymax/2).^2)/(2*radius^2));
S.fp0 = @(t) ones(size(t));
[S.fnvu, params] = nvu('REGDISABLE', false, 'CSTAR', cstar, 'GAMMA', gamma, 'FMET', fmet);

% Compute the new equilibrium
ua{i} = S.computeeq(u0, 1);
S.pressureflow(ua{i}, 1);
p{i} = S.tree.p;
q{i} = S.tree.q;
disp(i)
end

%% Generate a batch of graphs
close all
dirname = 'alphart2';
for i = 1:nr;
radius = Radii(i);
fmet = @(t, x, y) 1 + 0.5*exp(-((x - xmax/2).^2 + (y - ymax/2).^2)/(2*radius^2));
u = ua{i};
surfopts = {'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'none'};
axesopts = {'xlim', [min(x), max(x)], 'ylim', [min(y), max(y)], 'fontsize', 16};

% First graph: metabolic rate
figure(1)
clf


surf(X, Y, fmet(0, X, Y), surfopts{:})
set(gca, 'zlim', [1, 1.8], axesopts{:})
view([45, 30])
camlight
lighting phong
print('-djpeg', sprintf('%s/meta%d.jpg', dirname, i)); 

% Second graph: flow
figure(2)
clf
surf(X, Y, q{i}(iq), surfopts{:})
view([45, 30])
set(gca, axesopts{:})
camlight

lighting phong
print('-djpeg', sprintf('%s/flow%d.jpg', dirname, i)); 

% Third graph: tissue CO2 concentration
figure(3)
clf
surf(X, Y, ua{i}(S.imask + params.iT), surfopts{:})
view([45, 30])
set(gca, axesopts{:}, 'zlim', [cstar - 0.1, cstar + 1])
camlight

lighting phong
print('-djpeg', sprintf('%s/co2%d.jpg', dirname, i)); 

% Fourth graph: smooth muscle activation
figure(4)
clf
surf(X, Y, ua{i}(S.imask + params.iSMC), surfopts{:})
set(gca, axesopts{:}, 'zlim', [0, 1])
view(45, 15)
camlight

lighting phong
print('-djpeg', sprintf('%s/smc%d.jpg', dirname, i)); 
end
%% Data extraction
idx = S.imask; 
istate = idx(1, :);
iq1 = iq(1, :);

% SMC activation
styles = {'k-', 'k:', 'k-.', 'k--'};
figure(1)
for i = 1:nr
    fmet = @(t, x, y) 1 + 0.5*exp(-((x - xmin).^2 + (y - ymin).^2)/(2*Radii(i)^2));
    ff = fmet(0, x, ymin);
    plot(x, ff, styles{i}, 'linewidth', 1)
    hold on
    
end
set(gca, 'fontsize', 16, 'xlim', [min(x), max(x)])
hold off
ylabel('metabolic activity')
print -deps metabolic.eps

figure(2)
for i = 1:nr
    f = ua{i}(params.iSMC + istate);
    plot(x, f, styles{i}, 'linewidth', 1)
    hold on
end
hold off
set(gca, 'fontsize', 16, 'ylim', [0 1], 'xlim', [min(x), max(x)])
xlabel('x')
ylabel('smooth muscle activation')
print -deps smc.eps

figure(3)
for i = 1:nr
    ct = ua{i}(params.iT + istate);
    plot(x, ct, styles{i}, 'linewidth', 1)
    hold on
end
hold off
set(gca, 'fontsize', 16, 'ylim', [3, 4.5], 'xlim', [min(x), max(x)])
xlabel('x')
ylabel('Tissue [CO_2]')
print -deps co2.eps

figure(4)
for i = 1:nr
    qq = q{i}(iq1);
    plot(x, qq, styles{i}, 'linewidth', 1);
    hold on
end
hold off
set(gca, 'fontsize', 16, 'xlim', [min(x), max(x)])
xlabel('x')
ylabel('Blood flow')
print -deps flow.eps






%%
figure(1)
clf
subplot(2,2,1)
surf(ua(params.iSMC + S.imask), 'edgecolor', 'none');
set(gca, 'zlim', [0 1])
title('SMC activation')
subplot(2,2,2)
surf(ua(params.iR + S.imask), 'edgecolor', 'none');
title('radius')
subplot(2,2,3)

surf(S.tree.q(iq), 'edgecolor', 'none');
title('q')
subplot(2,2,4)
surf(ua(params.iT + S.imask), 'edgecolor', 'none');
title('Tissue CO_2')
set(gca, 'zlim', [cstar, cstar + 1]);

%% Hacking together some decent images
figure(1)
% Neurological activity
X = S.x(iq);
Y = S.y(iq);
x = x(1, :);
y = y(:, 1).';
imagesc(x, y, fmet(0, X, Y), [1 1.5])
colormap gray

figure(2)
imagesc(x, y, S.tree.q(iq))
colormap gray






