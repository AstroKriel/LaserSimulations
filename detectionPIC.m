function [INDICE, EXTREMA, INDICE2, MINIMUM, INDICE3, SS] = detectionPIC(DATA, j)
    INDICE = zeros(250,1);
    EXTREMA = zeros(250,1);
    compteur = 1;

    for i = (j-1)*8000+2:j*8000
        Courant = DATA(i);
        Precedent = Courant - DATA(i-1);
        Suivant = DATA(i+1) - Courant;
        
        % Si la valeur précédente et la valeur suivante sont inférieurx à la
        % valeur courante, on peut considérer qu'on se trouve sur un
        % extremum
        if(Precedent > 0 && Suivant <0 && compteur < 250)
            EXTREMA(compteur) = Courant;
            INDICE(compteur) = i;
            compteur = compteur + 1;
        end

    end

    EXTREMA = EXTREMA(1:compteur-1);
    INDICE = INDICE(1:compteur-1);

    INDICE2 = zeros(250,1);
    MINIMUM = zeros(250,1);
    compteur = 1;

    for i = (j-1)*8000+2:j*8000
        Courant = DATA(i);
        Precedent = Courant - DATA(i-1);
        Suivant = DATA(i+1) - Courant;

        % Si la valeur précédente et la valeur suivante sont inférieurx à la
        % valeur courante, on peut considérer qu'on se trouve sur un
        % extremum
        if(Precedent < 0 && Suivant > 0 && compteur < 250)
            MINIMUM(compteur) = Courant;
            INDICE2(compteur) = i;
            compteur = compteur + 1;
        end

    end

    MINIMUM = MINIMUM(1:compteur-1);
    INDICE2 = INDICE2(1:compteur-1);

    INDICE3 = zeros(250,1);
    SS = zeros(250,1);
    compteur = 1;
    epsi=10^(-8);

    for i = (j-1)*8000+2:j*8000
        Courant = DATA(i);
        Precedent = abs(Courant) - abs(DATA(i-1));
        Suivant = abs(DATA(i+1)) - abs(Courant);

        % Si la valeur précédente et la valeur suivante sont inférieurx à la
        % valeur courante, on peut considérer qu'on se trouve sur un
        % extremum
        if(abs(Precedent) < epsi && abs(Suivant) < epsi && compteur < 250)
            SS(compteur) = 1;
            INDICE3(compteur) = i;
            compteur = compteur + 1;
        end

    end

    SS = SS(1:compteur-1);
    INDICE3 = INDICE3(1:compteur-1);
end