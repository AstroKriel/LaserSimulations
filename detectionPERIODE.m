function [Tab_PERIODE, Tab_nb] = detectionPERIODE(DATA)
    max_i  = 12;
    Tab_PERIODE = zeros(max_i, 1);
    Tab_nb = zeros(max_i, 1);

    % On d�coupe la simulation (pour un th�ta et un gamma donn�s) en max_i parties
    % (8000 points) et on analyse la p�riodicit� du signal (fr�quence en nombre
    % de pics)
    for j = 1:max_i
        % Trouve les pics (les MAX uniquement)
        [INDICE, EXTREMA, ~, MINIMUM, ~, SS] = detectionPIC(DATA, j);

        % Il se peut que le signal soit constant ce qui implique que le vecteur
        % d'EXTREMA soit vide; auquel cas il n'y a pas de pic et donc aucune
        % p�riode. Dans le cas contraire, on lance l'algorithme d'analyse
        if isempty(EXTREMA) || isempty(MINIMUM) || abs(abs(EXTREMA(1))-abs(MINIMUM(1))) < 10^(-6)
            if SS  ==  1
                Tab_nb(j) = -10;
                Tab_PERIODE(j) = -10;
            else
                Tab_nb(j) = 0;
                Tab_PERIODE(j) = 0;
            end
        else
            % On soustrait le vecteur d'EXTREMA par le 1er EXTREMA trouv�.
            % Ainsi, si cet extrema est rep�r� plusieurs fois, le vecteur
            % r�sultant pr�sentera des "0" (valeurs des faibles) au niveau de
            % chaque extrema identique
            DIFF = EXTREMA - EXTREMA(1);
            INDICE_PIC = zeros(250, 1);
            compteur = 1;
            nb_periode = zeros(length(DIFF), 1);

            % On recherche la position des pics (donc des 0) en comparant chaque
            % terme de notre vecteur (DIFF) � un seuil de tol�rence (ici 0.0005)
            % Si un extremum est trouv�, on sauvegarde son indice (du vecteur de
            % simulation)
            % On conserve aussi sa position dans le vecteur d'EXTREMA qui sera
            % utile pour d�terminer le nombre de pics interm�diaires
            for i = 1:length(DIFF)
                if(abs(DIFF(i)) < 0.05)
                    INDICE_PIC(compteur) = INDICE(i);
                    nb_periode(compteur) = i;
                    compteur = compteur + 1;
                end
            end

            % INDICE_PIC contient l'indice (dans le vecteur de simulation) de
            % chaque extremum identique au premier trouv� (EXTREMA(1))
            INDICE_PIC = INDICE_PIC(1:compteur-1);
            % nb_periode contient l'emplacement (dans le vecteur EXTREMA) des
            % extremums identiques � l'EXTREMA(1) (cf. vecteur INDICE_PIC)
            nb_periode = nb_periode(1:compteur-1);

            % On �tudie ici la p�riodicit� (ou non) du signal
            INDICE_verif = zeros(compteur-1, 1);
            nb_periode_verif = zeros(compteur-1, 1);

            % Si le signal est p�riodique, l'�cart d'indice entre deux extrema
            % identiques doit �tre � peu pres constant. C'est ce qu'on v�rifie
            % ici par le biais du vecteur INDICE_verif
            % De plus, si le signal est p�riodique et que chaque p�riode
            % comprend plusieurs pics diff�rents, on doit pouvoir retrouver le
            % m�me nombre de pic dans chaque p�riode; c'est ce qu'on v�rifie
            % avec le vecteur nb_periode_verif
            for i = 2:length(INDICE_verif)
                INDICE_verif(i-1) = (INDICE_PIC(i) - INDICE_PIC(1))/(i-1);
                nb_periode_verif(i-1) = (nb_periode(i) - nb_periode(1))/(i-1);
            end

            nb_periode_verif = nb_periode_verif(1:length(INDICE_verif)-1);
            INDICE_verif = INDICE_verif(1:length(INDICE_verif)-1);

            % Si le signal est p�riodique, l'�cart d'indice entre deux extrema
            % identiques doit �tre plus ou moins constant. On analyse donc la
            % variance des �carts calcul�s pr�c�demment. Si elle est suffisement
            % faible, on peut consid�rer que la p�riode est la moyenne de chaque
            % �cart
            % Si la p�riodicit� est v�rifi�e, on stock le nombre de pics trouv�s
            % dans chaque p�riode dans le vecteur NB. Sinon, on consid�re qu'il
            % n'y a pas de p�riode ni de pics semblables (la plupart du temps
            % dans les zones de chaos)
            if round(var(INDICE_verif))  <  30
                PERIODE = abs(mean(INDICE_verif));
                NB = mean(nb_periode_verif);
            else
                PERIODE = 0;
                NB = 0;
            end

            % Il se peut que notre analyse trouve des p�riodes comprenant un
            % nombre de pics tr�s grands. On se limite � 5pics/p�riode pour se
            % consid�rer dans une zone de stabilit�. Si ce n'est pas le cas, on
            % dit que le signal est ap�riodique (PERIODE = 0)
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
