function res = descent(cx, alpha, t)

myf = @(x) x.^2;
mydf = @(x) 2*x;
alldata = [];
while abs(mydf(cx)) > 0.001
    del_x = -sign(mydf(cx));
    [xn, fn, fcall] = backtrack(cx, del_x, myf(cx), myf, mydf(cx), alpha, t, 0.001);
    alldata = [alldata; [xn, fn, fcall]];
    cx = xn;
end;

res = alldata;
    