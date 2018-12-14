[tab_periode,~] = detectionPERIODE(data);

if min(tab_periode) == -10 % stable
    C(index) = 0;
    C2(index) = 0;
    C3(index) = 0;
    is_stable(index) = true;
else
    if min(tab_periode) == 0 % apériodique
        C(index) = chaosBD(P,F); % La bande passante de chaos
        
        % C3(index) = eff_chaosBD(P,F);
        if C(index) < 10e6 % pour régler le problème des
            transitoires infinis de la restabilisation
            C(index) = 0;
            C2(index) = 0;
            C3(index) = 0;
            is_stable(index) = true;
        else
            coupe = find(F>2.3e10,1,'first');
            C2(index) = chaosBD(P(1:coupe),F(1:coupe)); % La
            bande passante de chaos
            is_chaos(index) = true;
        end
    else % périodique
        unevaleur = 1/(mean(tab_periode)*tau_p);
        C(index)=unevaleur;
        C3(index) = unevaleur;
        C2(index) = min(unevaleur,23e9);
        is_ECM(index)=true;
    end
end
