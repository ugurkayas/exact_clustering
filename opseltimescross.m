function y = opseltimescross(a,b)
y.val = a.val*b.val;
y.list = bagcrossjoin(b.list,a.list);
