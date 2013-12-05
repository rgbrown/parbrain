%% Autoregulation function design
% We will first consider the range of tissue concentrations that can be
% produced by our model

% First create a tree - maybe 13 levels is adequate
% 1/0.6 is good for RMULT

S = TreeSimulator('N', 13, 'RMULT', sqrt(2), 'NS', 8);
[fnvu, params] = nvu('REGDISABLE', true);
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
S.fnvu = nvu('REGDISABLE', true, 'CSTAR', cstar);
u0 = S.computeeq(u0);
iq = reshape(S.tree.i_spatial, sqrt(S.nBlocks), []);
X = S.x(iq);
Y = S.y(iq);
x = X(1, :);
y = Y(:, 1).';

%% Look at the range of q values that can be achieved by relaxing arterioles

% Look at increase in q by 

ua = u0;
for na = 1:5:100;
ua(params.nVars * (0:na-1) + params.iSMC) = 0;
[u0, success] = S.computeeq(ua);
S.pressureflow(u0, 1);
qa(na) = S.tree.q(1);
ct(na) = u0(params.iT);
disp(qa(na))
end

%% Now let's try with regulation on
% An important parameter to consider is GAMMA, this determines how steep
% the response is
[S.fnvu, params] = nvu('REGDISABLE', false, 'CSTAR', cstar, 'GAMMA', 100);

%% First, we'll look at changes in pressure - how well is it regulated
phigh = 1.5;
plow = 0.8;
uhigh = S.computeeq(u0, phigh);
S.pressureflow(uhigh, phigh);
qhigh = S.tree.q(1);
ulow = S.computeeq(u0, plow);
S.pressureflow(ulow, plow);
qlow = S.tree.q(1);



%% Now, we'll look at changes with neurological activity levels
% Create a metabolic function that starts in one corner, and we're just
% going to plot results for a cross section along one edge

Radii = [5e-4, 1e-3, 2e-3, 3e-3];
nr = numel(Radii);
for i = 1:nr

radius = Radii(i);
xmin = x(1); ymin = y(1);
fmet = @(t, x, y) 1 + 0.5*exp(-((x - xmin).^2 + (y - ymin).^2)/(2*radius^2));
[S.fnvu, params] = nvu('REGDISABLE', false, 'CSTAR', 3.5, 'GAMMA', 100, 'FMET', fmet);

% Compute the new equilibrium
ua{i} = S.computeeq(u0, 1);
S.pressureflow(ua{i}, 1);
p{i} = S.tree.p;
q{i} = S.tree.q;
end
%% Data extraction
idx = S.imask; 
istate = idx(1, :);
iq1 = iq(1, :);

% SMC activation
styles = {'k-', 'k:', 'k-.', 'k--'};
subplot(2,2,1);
for i = 1:nr
    fmet = @(t, x, y) 1 + 0.5*exp(-((x - xmin).^2 + (y - ymin).^2)/(2*Radii(i)^2));
    ff = fmet(0, x, ymin);
    plot(x, ff, styles{i}, 'linewidth', 2)
    hold on
end
set(gca, 'fontsize', 16)
hold off


subplot(2,2,2);
for i = 1:nr
    f = ua{i}(params.iSMC + istate);
    plot(x, f, styles{i}, 'linewidth', 2)
    hold on
end
hold off
set(gca, 'fontsize', 16, 'ylim', [0 1])
xlabel('x')
ylabel('f')

subplot(2,2,4);
for i = 1:nr
    ct = ua{i}(params.iT + istate);
    plot(x, ct, styles{i}, 'linewidth', 2)
    hold on
end
hold off
set(gca, 'fontsize', 16, 'ylim', [3, 4.5])
xlabel('x')
ylabel('c_t')

subplot(2,2,3)
for i = 1:nr
    qq = q{i}(iq1);
    plot(x, qq, styles{i}, 'linewidth', 2);
    hold on
end
hold off
set(gca, 'fontsize', 16)
xlabel('x')
ylabel('q')






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






