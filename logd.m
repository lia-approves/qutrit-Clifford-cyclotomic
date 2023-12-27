function logd = logd(n,b)
    global d;
    if nargin < 2
        b = d;
    end
    logd = round(log(n)/log(b),b);
end