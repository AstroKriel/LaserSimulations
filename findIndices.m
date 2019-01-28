function [itter_start, itter_end_1, itter_end_2] = findIndices(hor, ver)
%{
    Authors:    Neco Kriel, and
                Guillaume Bouchez (2019)

    Purpose:

    Input:
    - hor
    - ver

    Output:
    - itter_start
    - itter_end_2
%}

    points_found = false;
    itter_start = 1;
    while (itter_start < length(hor) && ~points_found)
        cur_hor = hor(itter_start);
        cur_ver = ver(itter_start);
        nex_hor = hor(itter_start + 1);
        nex_ver = ver(itter_start + 1);
        if (cur_hor <= cur_ver && nex_hor >= nex_ver && ~points_found)
            % hor grad > 0 & overlapping -> start period found
            itter_end_1 = itter_start + 1;
            while (itter_end_1 < length(hor) && ~points_found)
                cur_hor = hor(itter_end_1);
                cur_ver = ver(itter_end_1);
                nex_hor = hor(itter_end_1 + 1);
                nex_ver = ver(itter_end_1 + 1);
                if (cur_hor >= cur_ver && nex_hor <= nex_ver && ~points_found)
                    % hor grad < 0 & overlapping -> end period found
                    itter_end_2 = itter_end_1 + 1;
                    while (itter_end_2 < length(hor) && ~points_found)
                        cur_hor = hor(itter_end_2);
                        cur_ver = ver(itter_end_2);
                        nex_hor = hor(itter_end_2 + 1);
                        nex_ver = ver(itter_end_2 + 1);
                        if (cur_hor <= cur_ver && nex_hor >= nex_ver && ~points_found)
                            % hor grad > 0 & overlapping -> end period found
                            points_found = true;
                            break;
                        end

                        % not overlapping
                        itter_end_2 = itter_end_2 + 1;
                    end
                end

                % not overlapping
                itter_end_1 = itter_end_1 + 1;
            end
        end

        % not overlapping
        itter_start = itter_start + 1;
    end
end