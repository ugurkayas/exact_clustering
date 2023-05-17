function z = bagcrossjoin(x,y)
Nx = length(x);
Ny = length(y);
k = 1;
z = {};
for i = 1:Nx
    for j = 1:Ny
        z{k} = bagjoin(x{i},y{j});
        k = k + 1;
    end
end
