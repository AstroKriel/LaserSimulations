function [Chaos_Bandwidth ] = eff_chaosBD(P, F)
    % Compute effective chaos bandwidth (see Lin, Chao, Wu 2012)
    Tot_E = sum(P);
    Val = 0;
    Psort = sort(P,'descend');
    ifreq = 1; % index on the vector of frequencies
    
    while Val < 0.8*Tot_E
        Val = Val + Psort(ifreq);
        ifreq = ifreq + 1;
    end
    
    Chaos_Bandwidth= F(ifreq);
end