function [S, fode, u0] = treesim(varargin)
P = parseinputs(varargin{:});

% Construct the tree
S = Htree('N', P.N, 'RMULT', P.RMULT, ...
    'LRR', P.LRR, 'RMIN', P.RMIN);

nA = S.n - S.m;
iR   = 4*(0:nA-1) + 1; iSMC = 4*(0:nA-1) + 2;
iB   = 4*(0:nA-1) + 3; iT   = 4*(0:nA-1) + 4;

% Compute dimensionless parameter groups
E0 = P.EPASSIVE;
P0 = P.PROOT;
pcap = P.PCAP / P0;
R0 = P.RMIN;
L0 = P.LRR * P.RMIN;
Q0 = pi * R0^4 * P0 / (8 * P.MU * L0);
H0 = P.RMIN * P.HRR;
M0 = P.MNORMAL;
C0 = P.CIN;
regdisable = P.REGDISABLE;
icp = P.ICP / P.P0;
% Vessel wall parameters
a1 = E0 * P.T0 * P.RMIN / (P.ETA * R0);
a2 = P0 * P.RMIN * P.T0 ./ (P.ETA * H0);
a3 = P.RMIN / R0;
a4 = 1 - P.RSCALE;
a5 = P.EACTIVE/P.EPASSIVE - 1;
% CO2 stuff
Lmin = L0;
Acap = (2*Lmin)^3 / (P.RCAP / 2 + 1 / P.AVR);
Vt = Acap / P.AVR;
Vb = P.RCAP / 2 * Acap;
d1 = P.T0 * Q0 / Vb;
d2 = P.T0 * P.PW * Acap / Vb;
g1 = P.PW * P.T0 * Acap / Vt;
g2 = M0 * P.T0 / C0;
b1 = P.T0 / P.TAUM;
gamma = P.GAMMA;
cstar = P.CSTAR;
if isempty(P.FP0)
    fp0 = @p0;
else
    fp0 = P.FP0;
end
if isempty(P.FMET)
    fmet = @met;
else
    fmet = P.FMET;
end
fode = @evaluate;

% Set the conductance of the tree
S.setconductance((S.R / R0).^4 ./ (S.L / L0));
u0 = ones(4*nA, 1);
x = S.X(1:nA);
y = S.Y(1:nA);

    function du = evaluate(t, u)
        % Solve first
        S.setconductance(u(iR).^4, 1:nA);
        S.solve(fp0(t), pcap);
        % Compute transmural pressure
        pt = 0.5 * (S.p(S.i_parent(1:nA)) + pcap) - icp;
        du = rhs(t, u, x, y, S.q(1:nA), pt);
        
    end

    function du = rhs(t, u, x, y, q, pt)
        r  = u(iR); f  = u(iSMC);
        cb = u(iB); ct = u(iT);
        
        r0 = a3*(1 - a4*f);
        du = zeros(size(u));
        du(iR) = -a1*(1 + a5*f).*(r./r0 - 1) + a2*r.*pt;
        if ~regdisable
            du(iSMC) = -b1*(f - 1./(1 + exp(gamma*(ct - cstar))));
        else
            du(iSMC) = 0;
        end
        du(iB) = d1*q.*(1 - cb) + d2*(ct - cb);
        du(iT) = -g1*(ct - cb) + g2*fmet(t, x, y);
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
p.addParamValue('ICP', 3000'); % Pa
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
p.addParamValue('FMET', []);

p.parse(varargin{:});
params = p.Results;
end

function p = p0(t)
p = ones(size(t));

    %p = 1 + 0.5 ./ (1 + exp(-(t - 500)));
end

function m = met(t, x, y)
m = 1 + exp(-(x.^2 + y.^2) / ((200e-5)^2)) * exp(-((t - 100) / 20).^2);
end