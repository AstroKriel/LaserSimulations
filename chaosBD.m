function [Chaos_Bandwidth ] = chaosBD(P,F)
    Tot_E = sum(P);
    Val = 0;
    ifreq = 1; % index on the vector of frequencies
    
    while Val < 0.8*Tot_E
        Val = Val + P(ifreq);
        ifreq = ifreq + 1;
    end
    
    Chaos_Bandwidth= F(ifreq);
end