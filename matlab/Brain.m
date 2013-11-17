classdef Brain < handle
    properties
        % These properties are the workspace that is required to evaluate
        % the RHS and Jacobian of the ODE system
        regdisable
        iR
        iSMC
        iB
        iT
        iG
        A
        At
        m
        n
        G
        level
        b
        p
        q
        N
        l
        L0
        R0
        E0
        M0
        C0
        Lstar
        Rstar
        hstar
        Q0
        P0
        pcap
        a1
        a2
        a3
        a4
        a5
        b1
        d1
        d2
        g1
        g2
        gamma
        cstar  
        fp0
    end
   
    methods
        
        function S = Brain(varargin)
        %BRAIN constructs the data structure, and initialises the
        %nondimensional parameter groups, etc.
        
            P = parseinputs(varargin{:});
            S = inittree(S, P);
            % Set unscaled conductances
            S = setconductance(S, P, 0);
            S = solve(S, P.PROOT, P.PCAP);
            S = setconductance(S, P, 1);
            S = setparameters(S, P);
            
        end
        function dy = evaluate(S, t, y)
            nA = S.n - S.m;
            dy = zeros(size(y));
            r = y(S.iR);
            f = y(S.iSMC);
            cb = y(S.iB);
            ct = y(S.iT);
            % Update conductances of autoregulating vessels 
            S.G(S.iG) = r.^4 ./ S.l;
            S.solve(S.fp0(t), S.pcap);
            e = 1 + S.a5 .* f;
            r0 = S.a3 .* (1 - S.a4 .* f);
            pt = repmat(S.p(1:((S.m + 1)/2))', 2, 1) + S.pcap; pt = 0.5*pt(:);
            
            %pt = S.p(1:nA) - S.pcap;
            dy(S.iR)   = -S.a1 .* e .* (r ./ r0 - 1) + S.a2 .* r .* pt;
            if ~S.regdisable || (t < 0)
            dy(S.iSMC) = -S.b1 .* (f - 1 ./ (1 + exp(S.gamma .* (ct - S.cstar))));
            else
                dy(S.iSMC) = 0;
            end
            dy(S.iB)   = S.d1 .* S.q(1:nA) .* (1 - cb) + S.d2 .* (ct - cb);
            dy(S.iT)   = -S.g1 .* (ct - cb) + S.g2;
            %disp([pt(1), e(1), r(1), r0(1)]); 
        end
        
       
           
        function S = setparameters(S, P)
            nA = S.n - S.m;
            iA = 1:nA;
            S.iR = 4*(0:nA-1) + 1;
            S.iSMC = 4*(0:nA-1) + 2;
            S.iB = 4*(0:nA-1) + 3;
            S.iT = 4*(0:nA-1) + 4;
            S.iG = (S.n + 1) * (0:nA-1) + 1;
            S.E0 = P.EPASSIVE;
            S.P0 = P.PROOT;
            S.pcap = P.PCAP / S.P0;
            S.R0 = P.RMIN;
            S.L0 = P.LRR * P.RMIN;
            S.Q0 = pi * S.R0^4 * S.P0 / (8 * P.MU * S.L0);
            S.l = S.Lstar(iA) / S.L0;
            S.M0 = P.MNORMAL;
            S.C0 = P.CIN;
            S.regdisable = P.REGDISABLE;
            % Vessel wall parameters
            S.a1 = S.E0 * P.T0 * S.Rstar(iA) / (P.ETA * S.R0);
            S.a2 = S.P0 * S.Rstar(iA) * P.T0 ./ (P.ETA * S.hstar(iA));
            S.a3 = S.Rstar(iA) / S.R0;
            S.a4 = repmat(1 - P.RSCALE, nA, 1);
            S.a5 = repmat(P.EACTIVE/P.EPASSIVE - 1, nA, 1);
            % CO2 stuff
            Lmin = min(S.Lstar);
            Acap = (2*Lmin)^3 / (P.RCAP / 2 + 1 / P.AVR);
            Vt = Acap / P.AVR;
            Vb = P.RCAP / 2 * Acap;
            S.d1 = repmat(P.T0 * S.Q0 / Vb, nA, 1);
            S.d2 = repmat(P.T0 * P.PW * Acap / Vb, nA, 1);
            S.g1 = repmat(P.PW * P.T0 * Acap / Vt, nA, 1);
            S.g2 = repmat(S.M0 * P.T0 / S.C0, nA, 1);
            S.b1 = repmat(P.T0 / P.TAUM, nA, 1);
            S.gamma = repmat(P.GAMMA, nA, 1);
            S.cstar = repmat(P.CSTAR, nA, 1);
            if isempty(P.FP0)
                S.fp0 = @p0;
            else
                S.fp0 = P.FP0;
            end
            
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
        
        function S = inittree(S, P)
            S.A     = adjacency(P.N);
            S.At    = S.A.';
            [S.m, S.n]= size(S.A);
            S.level = floor(log2((S.n:-1:1)));
            S.b     = zeros(S.n, 1);
            S.p     = zeros(S.m, 1);
            S.q     = zeros(S.n, 1);
            S.Rstar = P.RMIN * P.RMULT.^(P.N - S.level - 1); 
            S.hstar = P.HRR * S.Rstar;
            S.Lstar = P.LRR * P.RMIN * 2.^(floor((P.N - S.level - 1) / 2));
            S.Rstar = S.Rstar(:); S.Lstar = S.Lstar(:); S.hstar = S.hstar(:);
            % Set the conductances
        end
        
        function S = setconductance(S, P, isscaled)
            if isscaled
                r = S.Rstar / P.R0;
                l = S.Lstar / P.L0; %#ok<PROP>
                S.G = spdiags(r.^4 ./ l, 0, S.n, S.n); %#ok<PROP>
            else
                S.G = spdiags(pi * S.Rstar.^4 ./ (8 * P.MU * S.Lstar), 0, S.n, S.n);
            end
        end        
    end
end
function params = parseinputs(varargin)
p = inputParser();
% Tree geometry
p.addParamValue('N', 9, @(x) x >= 3 && mod(x, 2) == 1); % odd int > 2
p.addParamValue('RMIN', 10e-6); %m
p.addParamValue('R0', 10e-6);
p.addParamValue('L0', 200e-6);
p.addParamValue('LRR', 20); %nondim
p.addParamValue('HRR', 0.1); %nondim
p.addParamValue('PROOT', 8000); % Pa
p.addParamValue('PCAP', 4000); % Pa
p.addParamValue('MU', 3.5e-3); % Pa s
p.addParamValue('RSCALE', 0.6); % nondim
p.addParamValue('RMULT', sqrt(2));
% Fluid flow
p.addParamValue('P0', 8000); % Pa
% Arterial wall
p.addParamValue('E0', 66e3); % Pa
p.addParamValue('EPASSIVE', 66e3); %Pa
p.addParamValue('EACTIVE', 233e3); %Pa
p.addParamValue('ETA', 2.8e2); % Pa s
%p.addParamValue('RHO', 1100); % kg m^-3
% Smooth muscle
p.addParamValue('GAMMA', 5); % nondim
p.addParamValue('CSTAR', 3.5); % nondim
p.addParamValue('TAUM', 5); % sec
% CO2
p.addParamValue('PW', 5e-5); %m s^-1
p.addParamValue('AVR', 5e3); %m^2 m^-3
p.addParamValue('RCAP', 5e-6); %m
p.addParamValue('CIN', 2.2e-2); %mol m^-3
p.addParamValue('MNORMAL', 4e-3); % mol m^-3 s^-1
% Timescale
p.addParamValue('T0', 1); %s
p.addParamValue('REGDISABLE', 0);
p.addParamValue('FP0', []);

p.parse(varargin{:});
params = p.Results;
end

 function p = p0(t)
    p = 1 + 0.5 ./ (1 + exp(-(t - 500)));
 end