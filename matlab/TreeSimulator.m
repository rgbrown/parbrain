classdef TreeSimulator < handle
    properties
        tree
        fp0
        fnvu
        nnvu
        nBlocks
        
        % These ones require reconstruction of the tree
        N
        RMULT
        
        % Required for initialising the tree
        R0
        L0
        PROOT
        PCAP
        P0
        PICP
        
        % Simulation stuff
        u0
        JPattern
        NS
        pcap
        picp
        x
        y
    end
    methods
        function S = TreeSimulator(varargin)
            P = parseinputs(varargin{:});
            fields = fieldnames(P);
            for i = 1:numel(fields)
                S.(fields{i}) = P.(fields{i});
            end
            S = inittree(S);
            S.setconductance();
        end
        function S = inittree(S)
            S.tree = Htree('N', S.N, 'RMULT', S.RMULT);
            S.nBlocks = S.tree.n - S.tree.m;
        end
        function setconductance(S)
            S.tree.setconductance((S.tree.R/S.R0).^4./(S.tree.L/S.L0));
            S.pcap = S.PCAP / S.P0;
            S.picp = S.PICP / S.P0;
            S.x = S.tree.X(1:S.nBlocks);
            S.y = S.tree.Y(1:S.nBlocks);
        end
        function setjpattern(S)
            Jb = sparse(ones(S.nnvu * S.NS));
            foo = repmat({Jb},S.nBlocks/S.NS,1);
            S.JPattern = blkdiag(foo{:});
        end
        function idx = imask(S)
            idx = 0:S.nnvu:(S.nnvu*(S.nBlocks - 1));
            idx = idx(S.tree.i_spatial);
            idx = reshape(idx, sqrt(S.nBlocks), []);
        end
        function du = evaluate(S, t, u, varargin)
            if nargin == 4
                p0 = varargin{1};
            else
                p0 = S.fp0(t);
            end
            % Solve first
            S.pressureflow(u, p0);
            % Compute transmural pressure
            pt = 0.5*(S.tree.p(S.tree.i_parent(1:S.nBlocks)) + S.pcap) ...
                - S.picp;
            du = S.fnvu(t, u, S.x, S.y, S.tree.q(1:S.nBlocks), pt);
        end
        function pressureflow(S, u, p0)
            S.tree.setconductance(u(1:S.nnvu:end).^4, 1:S.nBlocks);
            S.tree.solve(p0, S.pcap);
            
        end
        function [u0, success] = computeeq(S, varargin)
            dT = 200;
            success = false;
            maxits = 20;
            u0 = ones(S.nnvu * S.nBlocks, 1);
            p0 = 1;
            if nargin >= 2
                u0 = varargin{1};
            end
            if nargin >= 3
                p0 = varargin{2};
            end
            f = @(t, u) S.evaluate(t, u, p0);
            if isempty(S.JPattern)
                S.setjpattern();
            end
            opts = odeset('JPattern', S.JPattern);
            for i = 1:maxits
                [~, U] = ode15s(f, [0 dT], u0, opts);
                u0 = U(end, :).';
                if norm(f(0, u0), inf) / norm(u0, inf)  < 1e-4
                    success = true;
                    
                    break
                end
            end
            if success == false
                warning('Equilibrium not found');
            end
        end
  
        
        
        
        function set.N(S, val)
            S.N = val;
            if ~isempty(S.tree) && ~isempty(S.L0) && ~isempty(S.R0)
                S.inittree();
                S.setconductance();
            end
        end
        function set.RMULT(S, val)
            S.RMULT = val;
            if ~isempty(S.tree) && ~isempty(S.R0) && ~isempty(S.L0)
                S.inittree();
                S.setconductance();
            end
        end
        function set.R0(S, val)
            S.R0 = val;
            if ~isempty(S.L0) && ~isempty(S.tree)
                S.setconductance();
            end
        end
        function set.L0(S, val)
            S.L0 = val;
            if ~isempty(S.R0) && ~isempty(S.tree)
                S.setconductance();
            end
        end
        function set.NS(S, val)
            S.NS = val;
            if ~isempty(S.nnvu)
                S.setjpattern();
            end
        end
        
    end
    
    
    
end

function params = parseinputs(varargin)
p = inputParser();
p.addParamValue('N', 9);
p.addParamValue('RMULT', sqrt(2));
p.addParamValue('R0', 10e-6);  %m
p.addParamValue('L0', 200e-6); %m
p.addParamValue('PROOT', 8000); %Pa
p.addParamValue('PCAP', 4000); %Pa
p.addParamValue('PICP', 3000); %Pa
p.addParamValue('P0', 8000); %Pa
p.addParamValue('fp0', []);

p.addParamValue('NS', 2);

p.addParamValue('fnvu', []);
p.addParamValue('nnvu', []);
p.parse(varargin{:});
params = p.Results;
end