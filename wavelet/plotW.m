function plotW(w, del, T)

K = length(w);

offsets = zeros(1, K);
offsets(1) = 0;

if nargin > 2
    if length(T) == 1
        T = T * ones(1, K);
    end
end

for i = 2:K
    
    offsets(i) = min(w{i-1}) - del  - max(w{i});
    w{i} = w{i} - max(w{i}) + min(w{i-1}) - del;
    
    %     sig{i} = sig{i} - max(sig{i} - sig{i-1}) - del;
    
end

color = 'black';

t = 1:length(w{1});
plot(t, w{1}, color)
if nargin > 2
    line(t([1 end]), T(1)*[1 1], 'color', 'red')
    line(t([1 end]), -T(1)*[1 1], 'color', 'red')
end

for i = 2:K
    t = 1:length(w{i});
    line(t, w{i}, 'color', color)
    
    if nargin > 2
        if i <= length(T)
            line(t([1 end]), T(i)*[1 1] + offsets(i), 'color', 'red')
            line(t([1 end]), -T(i)*[1 1] + offsets(i), 'color', 'red')
        end
    end
end

% xlim([t(1) t(end)])

% M = round(0.8*length(sig{1}));
% for i = 1:K
%     V = max(sig{i}(M:end)) + 0.2*del;
%     th = text(t(end), V, labels{i}, 'horizontalAlignment', 'right', 'verticalalignment', 'bottom');
% end

% box off
% set(gca, 'ytick', [])

set(gca, 'tickdir', 'out')

m1 = max(w{1}) + 0*del;
m2 = min(w{K}) - 0*del;
ylim([m2 m1] + (m1-m2)*[-0.02 0.02])
