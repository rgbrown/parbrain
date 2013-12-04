%% Flow variability with recruitment and scaling
% This script demonstrates how the ability of a blood vessel to request
% extra blood by vasodilation is highly dependent on the morphology of the
% tree, and also and the amount of localised blood requested.
%
% We perform three experiments, the first shows the pressure distribution
% throughout the tree as a function of the scaling ratio of the vessels.
%
% For the sake of simplicity, we only look at scaling properties. Vessel
% lengths could be adjusted to give the same overall flow/resistance
% properties in each case.
%
% We consider a symmetric binary tree of $N$ levels, where the levels are
% numbered from the top down from $0$ to $N-1$. The radius $r_k$ of vessels
% at level $k$ is given by
% 
% $$r_k = \alpha^k r_0$$
% 
% And we assume that the length is a constant multiple of the radius, 
% 
% $$l_k = h_r r_k$$
%
% Up to a constant scale factor, then, the conductance of a vessel is 
% given by
%
% $$ g_k = r_k^{-3}$$
%

%% Preliminaries
% Initially, we construct an Htree object, and set the range of $\alpha$
% values to test.
N = 11;
H = Htree('N', N);
Alpha = linspace(1.1, 2, 10);
Rmin = 10e-6;
Hrr = 20;
idx = H.m - 2.^(0:N-2) + 1;

%%
% In this experiment we observe the relative increase in flow divided by
% the relative increase in radius, and see how these values can change
% markedly both with changes in alpha and the number of vessels recruited
dr = 1e-3;
sensitivity = zeros(numel(Alpha), N); 
P = zeros(N+1, numel(Alpha)); 
for i = 1:numel(Alpha)
    for j = 1:N
        % compute baseline pressure and flow
        [R0, L0, p0, q0] = initialise_vessels(H, Alpha(i), Rmin, Hrr, 1, 0);
    
        % increase the radius of the first 2^k vessels
        R = R0;
        R(1:2^(j-1)) = (1 + dr) * R(1:2^(j-1));
        H.setconductance(R.^4 ./ L0);
        H.solve(1, 0);
        sensitivity(i, j) = abs((H.q(1) - q0(1))/q0(1)) / ...
        abs((R(1) - R0(1))/R0(1));
    end
    P(:, i) = [1; p0(idx); 0];
end
figure(1)
plot(Alpha, P, 'k')
xlabel('\alpha')
ylabel('Pressure distribution')
figure(2)
plot(2.^(0:N-1), sensitivity)
set(gca, 'xscale', 'log')
xlabel('area activated')
legend(cellfun(@(s) strcat('\alpha = ', num2str(s)), num2cell(Alpha), 'uniformoutput', 0))
