A = zeros(5);
A(1,2) = 1; A(2,3) = 1;
A(2,5) = 1; A(3,4) = 1;
A = 100*eye(5) + A + A';

reachSets = cell(5, 1);

for i = 1:5
    [reachSets{i}, ~] = my_reach(A, i, 1:i-1);
end

for i = 1:5
    fprintf('Vertex %d: reach(%d, {%s}) = {%s}\n', i, i, num2str(1:i-1), num2str(reachSets{i}));
end

function [R, visited] = my_reach(A, v, S, R, visited)
    if nargin == 3
        R = [];
        visited(1:size(A, 2)) = false;
    end
    visited(v) = true;
    edges = find(abs(A(:, v)) > 0);
    if isempty(S)
        R = setdiff(edges, v);
        return;
    end
    for w = edges(:)'
        if ~visited(w)
            if ~ismember(S, w)
                R = [R w];
            else
                [R, visited] = my_reach(A, w, S, R, visited);
            end
        end
    end
end
