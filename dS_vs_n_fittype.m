function y = dS_vs_n_fittype( x, lambda, n_on_k, A, B )

y = n_on_k .* (A .* x.^lambda + B);

end
