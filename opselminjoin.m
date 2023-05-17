function y = opselminjoin(a,b)
if (a.val < b.val)
    y = a;
elseif (a.val > b.val)
    y = b;
else
    y.val = a.val;
    y.list = bagjoin(a.list,b.list);
end
