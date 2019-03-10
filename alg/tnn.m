function y = tnn(X, r)

sigs = svd(X);

y = sum(sigs(r+1:end));

end