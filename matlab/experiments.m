%% Brain code test script

% Influence of radius multiple on overall tree diameter

mult = linspace(1, 2, 20);
for i = 1:numel(mult)
    S = Brain('N', 13, 'RMULT', mult(i));
    S.solve(1, 0.5);
    q(i) = S.q(end) * S.Q0;
end
plot(mult, q)

%% Test maximum and minimum concentrations for different values of rmult
mult = sqrt(2)


for i = 1:numel(mult)
    % Construct tree and Jacobian
    S = Brain('N', 13, 'CSTAR', 2, 'GAMMA', 5, 'FP0', @(t) ones(size(t)), ...
        'REGDISABLE', true, 'RMULT', mult(i));
    nvars = (S.n - S.m) * 4;
    f = @(t, x) S.evaluate(t, x);
    k = 2;
    clear foo
    Jb = sparse(ones(4*k));
    for i = (1:(S.n - S.m) / k)
        foo{i} = Jb;
    end
    J = blkdiag(foo{:});
    opts = odeset('JPattern', J);
    % Run to steady state
    y0 = ones(4*(S.n - S.m), 1);
    y0(S.iSMC) = 1;
    dT = 50;
    while norm(f(0, y0), inf) > 1e-6
        [~, X] = ode15s(f, [0 dT], y0, opts);
        y0 = X(end, :).';
       
    end
    disp('y0 found');

    [T, Y] = ode15s(f, [0 200], y0, opts);
    S.compute_flow(T(end), Y(end, :).');
    Qmin(i) = S.q(1);
    cmin(i) = max(Y(end, S.iT));
    
    
    y0 = ones(4*(S.n - S.m), 1);
    y0(S.iSMC) = 0;
    dT = 200;
    while norm(f(0, y0), inf) > 1e-6
        [~, X] = ode15s(f, [0 dT], y0, opts);
        y0 = X(end, :).';
        
    end
    disp('y0 found');
    
    [T, Y] = ode15s(f, [0 200], y0, opts);
    S.compute_flow(T(end), Y(end, :).');
    cmax(i) = min(Y(end, S.iT));
    Qmax(i) = S.q(1);
    disp([cmin(i), cmax(i)]);
    disp([Qmin(i), Qmax(i)]);
    
    
    
    
    
    
    
    
end


%%
mult = 1.6
S = Brain('N', 13, 'CSTAR', 3, 'GAMMA', 5, 'FP0', @(t) ones(size(t)), ...
    'REGDISABLE', true, 'RMULT', 1.6);





