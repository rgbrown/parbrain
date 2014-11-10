function [R_base, L_base, p_base, q_base] = initialise_vessels(H, alpha, Rmin, Hrr, p0, pcap)
    N = H.N;
    r = Rmin * alpha.^(0:N-1);
    R_base = r(N - H.level);
    L_base = Hrr * R_base;
    H.R = R_base;
    H.L = L_base;
    H.setconductance(H.R.^4 ./ H.L);
    H.solve(p0, pcap);
    p_base = H.p;
    q_base = H.q;
end