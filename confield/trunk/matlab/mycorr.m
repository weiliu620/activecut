function RHO = mycorr(X, Y)
% Usage: mycorr(X, Y). my version fo matlab's corr function.

if size(X, 1) ~= size(Y,1)
    error(mycorr:dim, 'X and Y must have same row');
end;
[N, P1] = size(X);
P2 = size(Y, 2);
RHO = zeros(P1, P2);
X = X - repmat(mean(X, 1), N, 1);
Y = Y - repmat(mean(Y, 1), N, 1);
for p1 = 1:P1
    for p2 = 1:P2
        RHO(p1,p2) = dot(X(:,p1), Y(:,p2))/((N-1)*std(X(:,p1))* std( Y(:,p2)));
    end;
end;