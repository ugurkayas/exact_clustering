function [y,c]= Exact_dp_means(x,lambda)

x=sort(x);                                         % sorting the data
N= length(x);                                      % length of the data
mu=zeros(N,N);                                     % global mean matris for each possible mean from i to j
eij = zeros(N,N);                                  % objective value matris for each possible mean from i to j
mu(1,1)=x(1);                                      % initialization of the mean matris
eij(1,1)=lambda;                                   % initialization of the objective value matris 
                                                                        
for j=(2:N)                                        % objective value matris calculation
    for i=(1:j)
        mu(i,j) = (1/(j-i+1))*((j-i)*mu(i,j-1)+x(j));
        if i==j
             eij(i,j) = lambda;
        else
             eij(i,j)= eij(i,j-1)  + ((j-i)/(j-i+1))*(x(j)-mu(i,j-1)).^2;
        end
    end
end


                                                   % min-plus semiring with tuples
opplus = @opselminjoin;                                                                 
idplus.val = inf;
idplus.list = {};
optimes = @opselpluscross;
idtimes.val = 0;
idtimes.list = {{}};

E = cell(N+1,1);                                   % Bellman's recursion
for i = 1:N+1
     E{i,1} = idplus;
end
E{1,1} = idtimes;
for j = 1:N
    E{j+1} = idplus;
    for i = 1:j
        E{j+1} = opplus(E{j+1},optimes(E{i},opselinit(eij(i,j),[i,j])));
    end
end

clusters = E{N+1,1}.list;                         % extracting the final clustering setting
for k = 1:length(clusters)
  c=sprintf('(%d,%d)',clusters{k}{:});            
end
c=append('[',c,']');                              % the final clustering setting                                       
y=E{N+1,1}.val;                                   % objective value

end
