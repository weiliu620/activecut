function testouter()

A = [1 2 2 3]';

[cell, mat] = outer(@kernel2, A, A)

function res = kernel(x,y)
    res = (x == y);
return

kernel2 = @(x,y) (x==y);
