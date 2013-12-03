%% Ratio Test

Rmin = 10e-6;
Rmax = 200e-6;
alpha = 1.6;

% We can compute N automatically, or just set it
% N = floor(log(Rmax / Rmin) / log(alpha));
N = 7;

eta = 2.8e2;
H = Htree('N', N, 'RMULT', alpha);
r = H.RMIN * alpha.^(0:H.N-1);
H.R = r(N - H.level);
H.L = 20 * H.R;
%%
g = H.R.^4 * pi ./ (8*eta*H.L);
H.setconductance(g);
p0 = 4000;
H.solve(p0, 0);

idx = H.m - 2.^(0:N-2) + 1;
r = H.R(H.n - 2.^(0:N-1) + 1);

plot(1e6*r, [p0; H.p(idx)], 'k', 'linewidth', 2)
set(gca, 'xdir', 'reverse')
xlabel('radius (\mu m)')
ylabel('pressure at head of vessel')
conductance = H.q(end) / p0;




% We want to create trees with different scaling, and see how the overall
% resistance varies
