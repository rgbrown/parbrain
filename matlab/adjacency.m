function A = adjacency(N)
%ADJACENCY   Compute adjacency matrix for N level tree
%    A = ADJACENCY(N)
%        N is the number of levels of branches in the tree
m = 2^(N-1) - 1;
n = 2^N - 1;
A = sparse(1:m, 2*(1:m)-1, 1, m, n, m) + ...
    sparse(1:m, 2*(1:m), 1, m, n, m) + ...
    sparse(1:m, (1:m) + (n - m), -1, m, n, m);
end
