function [protrusion_lifetime, protrusion_len, protrusion_P_avg] = GetProtrusion(Pangle_all_2cell,P_all_2cell, ang_thres)
% Protrusion %
dang_thres = ang_thres/180*pi();               % threshold of change in angle
protrusion_counter = 1;
idx_i = 1;
idx_f = [];

Pangle_all_2cell_tmp = Pangle_all_2cell;

% Find start and end time of protrusions
while length(Pangle_all_2cell_tmp) >= 1
    ang_current = Pangle_all_2cell_tmp(1);
    idx = find(abs(angdiff(Pangle_all_2cell_tmp,ang_current*ones(size(Pangle_all_2cell_tmp)))) >= dang_thres,1);
    if isempty(idx)
        idx_f(protrusion_counter) = length(Pangle_all_2cell);
        break
    end
    idx_i(protrusion_counter + 1) = idx_i(protrusion_counter) + idx - 1;
    idx_f(protrusion_counter) = idx_i(protrusion_counter + 1) - 1;
    Pangle_all_2cell_tmp(1:(idx-1)) = [];
    protrusion_counter = protrusion_counter + 1;
end

[~, protrusion_P_avg, ~] = deal(zeros(1,protrusion_counter));
protrusion_lifetime = idx_f-idx_i;

for jj = 1:protrusion_counter
    protrusion_P_avg(jj) = mean(P_all_2cell(idx_i(jj):idx_f(jj)));
end
protrusion_len = protrusion_P_avg.*protrusion_lifetime;



end