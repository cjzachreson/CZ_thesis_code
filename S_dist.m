function [S] = S_dist(vars, minimum, maximum, n_bins)

range = maximum - minimum;

edges = [minimum:range/n_bins:maximum];

pop_dist = histogram(vars, edges, 'Visible', 'off');

counts = pop_dist.Values;

binWidth = range/n_bins;

nz = counts > 0;

px = counts(nz) / sum(counts(nz));

S = -sum(px.*log(px));

end

