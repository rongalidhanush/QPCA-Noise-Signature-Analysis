function bins = select_qpca_peaks(probability,K)
%SELECT_QPCA_PEAKS Select K separated local maxima from a QPE spectrum.
%
% The function prevents adjacent bins belonging to the same finite-
% resolution QPE peak from being counted as separate components.

p = probability(:);
N = numel(p);

local = false(N,1);

for i = 1:N
    left  = p(mod(i-2,N)+1);
    right = p(mod(i,N)+1);
    local(i) = (p(i) >= left) && (p(i) >= right);
end

candidates = find(local);
[~,ord] = sort(p(candidates),'descend');
candidates = candidates(ord);

selected = [];
minimumSeparation = max(1,round(N/64));

for i = 1:numel(candidates)
    c = candidates(i);
    if isempty(selected)
        selected(end+1) = c; %#ok<AGROW>
    else
        circularDistance = min(abs(c-selected),N-abs(c-selected));
        if all(circularDistance >= minimumSeparation)
            selected(end+1) = c; %#ok<AGROW>
        end
    end

    if numel(selected) == K
        break;
    end
end

% Fallback if the spectrum has fewer local maxima than requested.
if numel(selected) < K
    [~,allOrd] = sort(p,'descend');
    for i = 1:numel(allOrd)
        c = allOrd(i);
        if isempty(selected) || all(min(abs(c-selected),N-abs(c-selected)) >= minimumSeparation)
            selected(end+1) = c; %#ok<AGROW>
        end
        if numel(selected)==K
            break;
        end
    end
end

bins = selected(:);

end
