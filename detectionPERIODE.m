function [Tab_PERIODE, Tab_nb] = detectionPERIODE(DATA)
    max_i  = 12;
    Tab_PERIODE = zeros(max_i, 1);
    Tab_nb = zeros(max_i, 1);

    % On découpe la simulation (pour un théta et un gamma donnés) en max_i parties
    % (8000 points) et on analyse la périodicité du signal (fréquence en nombre
    % de pics)
    for j = 1:max_i
        % Trouve les pics (les MAX uniquement)
        [INDICE, EXTREMA, ~, MINIMUM, ~, SS] = detectionPIC(DATA, j);

        % Il se peut que le signal soit constant ce qui implique que le vecteur
        % d'EXTREMA soit vide; auquel cas il n'y a pas de pic et donc aucune
        % période. Dans le cas contraire, on lance l'algorithme d'analyse
        if isempty(EXTREMA) || isempty(MINIMUM) || abs(abs(EXTREMA(1))-abs(MINIMUM(1))) < 10^(-6)
            if SS  ==  1
                Tab_nb(j) = -10;
                Tab_PERIODE(j) = -10;
            else
                Tab_nb(j) = 0;
                Tab_PERIODE(j) = 0;
            end
        else
            % On soustrait le vecteur d'EXTREMA par le 1er EXTREMA trouvé.
            % Ainsi, si cet extrema est repéré plusieurs fois, le vecteur
            % résultant présentera des "0" (valeurs des faibles) au niveau de
            % chaque extrema identique
            DIFF = EXTREMA - EXTREMA(1);
            INDICE_PIC = zeros(250, 1);
            compteur = 1;
            nb_periode = zeros(length(DIFF), 1);

            % On recherche la position des pics (donc des 0) en comparant chaque
            % terme de notre vecteur (DIFF) à un seuil de tolérence (ici 0.0005)
            % Si un extremum est trouvé, on sauvegarde son indice (du vecteur de
            % simulation)
            % On conserve aussi sa position dans le vecteur d'EXTREMA qui sera
            % utile pour déterminer le nombre de pics intermédiaires
            for i = 1:length(DIFF)
                if(abs(DIFF(i)) < 0.05)
                    INDICE_PIC(compteur) = INDICE(i);
                    nb_periode(compteur) = i;
                    compteur = compteur + 1;
                end
            end

            % INDICE_PIC contient l'indice (dans le vecteur de simulation) de
            % chaque extremum identique au premier trouvé (EXTREMA(1))
            INDICE_PIC = INDICE_PIC(1:compteur-1);
            % nb_periode contient l'emplacement (dans le vecteur EXTREMA) des
            % extremums identiques à l'EXTREMA(1) (cf. vecteur INDICE_PIC)
            nb_periode = nb_periode(1:compteur-1);

            % On étudie ici la périodicité (ou non) du signal
            INDICE_verif = zeros(compteur-1, 1);
            nb_periode_verif = zeros(compteur-1, 1);

            % Si le signal est périodique, l'écart d'indice entre deux extrema
            % identiques doit être à peu pres constant. C'est ce qu'on vérifie
            % ici par le biais du vecteur INDICE_verif
            % De plus, si le signal est périodique et que chaque période
            % comprend plusieurs pics différents, on doit pouvoir retrouver le
            % même nombre de pic dans chaque période; c'est ce qu'on vérifie
            % avec le vecteur nb_periode_verif
            for i = 2:length(INDICE_verif)
                INDICE_verif(i-1) = (INDICE_PIC(i) - INDICE_PIC(1))/(i-1);
                nb_periode_verif(i-1) = (nb_periode(i) - nb_periode(1))/(i-1);
            end

            nb_periode_verif = nb_periode_verif(1:length(INDICE_verif)-1);
            INDICE_verif = INDICE_verif(1:length(INDICE_verif)-1);

            % Si le signal est périodique, l'écart d'indice entre deux extrema
            % identiques doit être plus ou moins constant. On analyse donc la
            % variance des écarts calculés précédemment. Si elle est suffisement
            % faible, on peut considérer que la période est la moyenne de chaque
            % écart
            % Si la périodicité est vérifiée, on stock le nombre de pics trouvés
            % dans chaque période dans le vecteur NB. Sinon, on considère qu'il
            % n'y a pas de période ni de pics semblables (la plupart du temps
            % dans les zones de chaos)
            if round(var(INDICE_verif))  <  30
                PERIODE = abs(mean(INDICE_verif));
                NB = mean(nb_periode_verif);
            else
                PERIODE = 0;
                NB = 0;
            end

            % Il se peut que notre analyse trouve des périodes comprenant un
            % nombre de pics très grands. On se limite à 5pics/période pour se
            % considérer dans une zone de stabilité. Si ce n'est pas le cas, on
            % dit que le signal est apériodique (PERIODE = 0)
            if NB  <  5
                Tab_nb(j) = NB;
                Tab_PERIODE(j) = PERIODE;
            else
                Tab_nb(j) = 0;
                Tab_PERIODE(j) = 0;
            end
        end
    end
end
