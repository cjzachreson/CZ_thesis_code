function [X_PB] = wrap_point(X, L)

if abs(X) > L/2
    X_PB = (L - abs(X)) * -sign(X);
    
else X_PB = X;
end
end