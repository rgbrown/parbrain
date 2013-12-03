%% Construct a Htree
N = 11;
H = Htree('N', N);

%alpha = 2^(1/3);
alpha = 1.3;
Rmin = 10e-6;
r = Rmin * alpha.^(0:N-1);
R0 = r(N - H.level);
H.R = R0;
H.L = 20 * H.R;

H.setconductance(H.R.^4 ./ H.L);

idx = H.m - 2.^(0:N-2) + 1;

% Compute baseline pressure/flow distribution
H.solve(1, 0);
qbase = H.q(1);
p0 = H.p;
q0 = H.q;

dr = 5e-1; % relative change
%% Single vessel perturbation experiment
a = 6;
% Now, perturb one of the radii
R = R0;
rbase = R(1);
rnew = (1 + dr) * rbase;
R([1, 1+ randperm(H.n-1, 2^a - 1)]) = rnew;
H.R = R;
H.setconductance(H.R.^4 ./ H.L);
H.solve(1, 0);
qnew = H.q(1);

sensitivity1 = abs((qnew - qbase)/qbase) / abs((rnew - rbase)/rbase)

%% All vessels perturbation experiment
R = R0;
rbase = R(1);
rnew = (1 + dr) * rbase;
R(H.level == N - 1) = rnew;
H.R = R;
H.setconductance(H.R.^4 ./ H.L);
H.solve(1, 0);
qnew = H.q(1);
sensitivityall = abs((qnew - qbase)/qbase) / abs((rnew - rbase)/rbase)


