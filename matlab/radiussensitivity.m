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
Alpha = linspace(0.1, 0.9, 5);
Rmin = 10e-6;
Hrr = 20;
idx = H.m - 2.^(0:N-2) + 1;

%%
% In this experiment we observe the relative increase in flow divided by
% the relative increase in radius, and see how these values can change
% markedly both with changes in alpha and the number of vessels recruited
dr = 0.25;
sensitivity = zeros(numel(Alpha), 2^(N-1)); 
P = zeros(N+1, numel(Alpha)); 
for i = 1:numel(Alpha)
    fprintf('%d ', i);
    for j = 1:2^(N-1)
        % compute baseline pressure and flow
        [R0, L0, p0, q0] = initialise_vessels(H, 1/Alpha(i), Rmin, Hrr, 1, 0);
    
        % increase the radius of the first 2^k vessels
        R = R0;
        R(1:j) = (1 + dr) * R(1:j);
        H.setconductance(R.^4 ./ L0);
        H.solve(1, 0);
        sensitivity(i, j) = abs((H.q(1) - q0(1))/q0(1));% / ...
        %abs((R(1) - R0(1))/R0(1));
    end
    P(:, i) = [1; p0(idx); 0];
end
fprintf('\n')
%%
% Set Alpha = linspace(0,1,50) and then run block above to do this

figure(1)
plot(Alpha, P, 'k', 'linewidth', 2)
set(gca, 'fontsize', 16)
set(gca, 'ylim', [0 1])
xlabel('\alpha')
ylabel('Normalised pressure')

% annotate
[~, idx] = min(abs(Alpha - 0.6));
text(Alpha(idx) - 0.07, P(end-1, idx), sprintf('p_{%d}', N-1), 'fontsize', 16)
[~, idx] = min(abs(Alpha - 0.85));
text(Alpha(idx) + 0.01, P(2, idx), 'p_1', 'fontsize', 16)
print -deps pressure.eps


%%
% Set Alpha = linspace(1.1, 1.9, 5) and run block above first
styles ={'k-', 'k:', 'k-.', 'k--'};
figure(2)
clf
for i = 1:size(sensitivity, 1)
    plot(1:2^(N-1), 100*sensitivity(i, :), styles{mod(i-1, 4)+1}, 'linewidth', 2)
hold on
end
    set(gca, 'fontsize', 16);
set(gca, 'xscale', 'log', 'xlim', [0, 1024])
xlabel('no. vessels activated')
ylabel('flow increase through single vessel (%)')
legend(cellfun(@(s) strcat('\alpha = ', num2str(s)), num2cell(Alpha), 'uniformoutput', 0))
print -deps recruitment.eps

