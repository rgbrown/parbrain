function [frhs, P] = nvu(varargin)
P = parseinputs(varargin{:});

% Vessel wall model parameters
E0 = P.EPASSIVE;
L0 = P.LRR * P.RMIN;
H0 = P.RMIN * P.HRR;
a1 = E0 * P.T0 * P.RMIN / (P.ETA * P.R0);
a2 = P.P0 * P.RMIN * P.T0 ./ (P.ETA * H0);
a3 = P.RMIN / P.R0;
a4 = 1 - P.RSCALE;
a5 = P.EACTIVE/P.EPASSIVE - 1;

% CO2 model
P.M0 = P.MNORMAL;
P.C0 = P.CIN;
Lmin = L0;
Acap = (2*Lmin)^3 / (P.RCAP / 2 + 1 / P.AVR);
Vt = Acap / P.AVR;
Vb = P.RCAP / 2 * Acap;
if isempty(P.FMET)
    P.FMET = @met;
end
fmet = P.FMET;
P.Q0 = pi * P.R0^4 * P.P0 / (8 * P.MU * L0);
d1 = P.T0 * P.Q0 / Vb;
d2 = P.T0 * P.PW * Acap / Vb;
g1 = P.PW * P.T0 * Acap / Vt;
g2 = P.M0 * P.T0 / P.C0;

% SMC model parameters
regdisable = P.REGDISABLE;
b1 = P.T0 / P.TAUM;
gamma = P.GAMMA;
cstar = P.CSTAR;

iR = 1;
iB = 2;
iT = 3;
iSMC = 4;
P.iR = iR; P.iB = iB; P.iT = iT; P.iSMC = iSMC;
P.nVars = 4;

frhs = @rhs;

    function du = rhs(t, u, x, y, q, pt)
        r  = u(iR:4:end);
        f  = u(iSMC:4:end);
        cb = u(iB:4:end);
        ct = u(iT:4:end);
        
        r0 = a3*(1 - a4*f);
        du = zeros(size(u));
        du(iR:4:end) = -a1*(1 + a5*f).*(r./r0 - 1) + a2*r.*pt;
        if ~regdisable
            du(iSMC:4:end) = -b1*(f - 1./(1 + exp(gamma*(ct - cstar))));
        else
            du(iSMC:4:end) = 0;
        end
        du(iB:4:end) = d1*q.*(1 - cb) + d2*(ct - cb);
        du(iT:4:end) = -g1*(ct - cb) + g2*fmet(t, x, y);
    end


end

function params = parseinputs(varargin)
p = inputParser();
% Tree geometry

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
p.addParamValue('FMET', []);

p.parse(varargin{:});
params = p.Results;

end

function m = met(t, x, y)
m = 1 + exp(-(x.^2 + y.^2) / ((200e-5)^2)) * ...
    exp(-((t - 100) / 50).^2);
end
