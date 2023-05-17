function y = opselpluscross(a,b)
y.val = a.val+b.val;
y.list = bagcrossjoin(a.list,b.list); % Grow start of bag
% y.list = bagcrossjoin(b.list,a.list); % Grow end of bag

