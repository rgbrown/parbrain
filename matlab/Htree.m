classdef Htree < handle
    properties
        A
        At
        R
        L
        G
        X
        Y
        p
        q
        g
        is_terminal
        i_spatial
        i_parent
        
        RMIN
        N
        b
        m
        n
        level
    end
    
    methods
        function S = Htree(varargin)
            params = parseinputs(varargin{:});
            S = inittree(S, params);
        end
        
        function S = inittree(S, params)
            S.N = params.N;
            S.A = adjacency(S.N);
            S.At = S.A.';
            [S.m, S.n]= size(S.A);
            S.level = floor(log2((S.n:-1:1)));
            S.b     = zeros(S.n, 1);
            S.p     = zeros(S.m, 1);
            S.q     = zeros(S.n, 1);
            S.RMIN = params.RMIN;
            S.R = params.RMIN * params.RMULT.^(params.N - S.level - 1); 
            S.L = params.LRR * params.RMIN * ...
                2.^(floor((params.N - S.level - 1) / 2));
            S.R = S.R(:);
            S.L = S.L(:);
            S.is_terminal = S.level == max(S.level);
            S = initspatial(S);
            S.G = speye(S.n);
            
            % Set the conductances
        end
        
        function S = solve(S, pin, pcap)
            S.b(end) = pin;
            S.b(1:(S.n - S.m)) = -pcap;
            AG = S.A*S.G;
            bh = AG * S.b;
            B = AG*S.At;
            clear AG;            
            % Our matrix is ordered such that it has a perfect Cholesky
            % factorisation, no column reordering required
            K = chol(B, 'lower');
            clear B;
            S.p = -K' \ (K \ bh);
            S.q = S.G*(S.At*S.p + S.b);
        end
        
        function S = setconductance(S, g, varargin)
            if nargin == 2
                idx = 1:S.n;
            else
                idx = varargin{1};
            end
            if islogical(idx)
                idx = find(idx);
            end
            S.G( (idx - 1) * (S.n + 1) + 1) = g;
        end
        
    end
    
end

function S = initspatial(S)
a = (1:(2^(S.N-1))) - 1;
x = zeros(size(a));
y = zeros(size(a));

for i = 0:S.N-2
    d = 2^(floor(i/2)) * (2*mod(a, 2) - 1);
    if ~mod(i, 2)
        x = x + d;
    else
        y = y + d;
    end
    a = bitshift(a, -1);
end

% work out spatial ordering
xmin = min(x);
ymin = min(y);
xs = (x - xmin) / 2;
ny = 2^floor((S.N-1)/2);
ys = (y - ymin) / 2;
idx = ys + ny * xs + 1;
S.i_spatial(idx) = 1:(S.n - S.m);
% Find the remainders by averaging x and y coordinates of children
A = S.A;
[m, n] = size(A);
A(A == -1) = -2;
A = [speye(n-m, n); A];

b = zeros(n, 1);
b(1:(n-m)) = x;
Lmin = min(S.L);
S.X = Lmin * (A \ b);
b(1:(n-m)) = y;
S.Y = Lmin * (A \ b);

% Now compute parents
S.i_parent = zeros(1, n);
[I, J] = find(S.At);
S.i_parent(I(1:3:end)) = J(3:3:end);
S.i_parent(I(2:3:end)) = J(3:3:end);
end

function params = parseinputs(varargin)
p = inputParser();
p.addParamValue('N', 9)
p.addParamValue('RMIN', 10e-6)
p.addParamValue('RMULT', sqrt(2))
p.addParamValue('LRR', 20)
p.parse(varargin{:});
params = p.Results;
end
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
