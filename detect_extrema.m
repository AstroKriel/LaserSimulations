function EXTREMA = detect_extrema(DATA, num_extrema, epsilon)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

params = [num_extrema epsilon length(DATA)];
[EXTREMA count] = fast_detect_extrema(DATA,params);
count = count+1;


% EXTREMA = zeros(num_extrema,1);
% 
% inc = 2;
% count = 1;
% 
% while (inc < length(DATA)) && (count <= num_extrema)
%     current = DATA(inc);
%     previous = sign(current-DATA(inc-1));
%     next = sign(DATA(inc+1)-current);
%     
%     if (previous ~= next)
%         if ~any(abs(EXTREMA-current) <= epsilon)
%             EXTREMA(count) = current;
%             count = count+1;
%         end
%     end
%     
%     inc = inc+1;
% end
% 
% if count == 1
%     EXTREMA(count) = DATA(inc);
%     count = count+1;
% end

EXTREMA = EXTREMA(1:count-1);

end

