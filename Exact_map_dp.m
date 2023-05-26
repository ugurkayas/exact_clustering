clear all;
clc;

x=[5.3,1,5.1,5.4,1.2,1.4,1.3,5,5.2,1.1];                     % sample data
x=sort(x);                                                   % sorting the data
N= length(x);                                                % length of the data

alpha0=0.05;                                                 % hyper-parameters
sigma0=(1)*std(x);
sigma_hat=(0.5)*std(x);
m0=mean(x);

for j = 1:N                                                  % objective value matris calculation
    for i = 1:j
        sigma_k=(1/sigma0^2+length(x(i:j))/sigma_hat^2)^(-1);
        mu_k=sigma_k*(m0/sigma0^2+sum(x(i:j))/sigma_hat^2);
        eij(i,j) = 0;
        for q=i:j
            eij(i,j)= eij(i,j) + Gaussian_nll(x(q),mu_k,sigma_k,sigma_hat);
        end
        eij(i,j) = eij(i,j) - log(alpha0) - gammaln(length(x(i:j)));
    end
end

                                                            % min-plus semiring with tuples
opplus = @opselminjoin;
idplus.val = inf;
idplus.list = {};
optimes = @opselpluscross;
idtimes.val = 0;
idtimes.list = {{}};

E = cell(N+1,1);                                            % Bellman's recursion
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

clusters = E{N+1,1}.list;                                    % extracting the final clustering setting
for k = 1:length(clusters)
    fprintf('[');
    fprintf('(%d,%d)',clusters{k}{:});
    fprintf(']');
end

fprintf('\n\n');                                            % the final clustering setting
E{N+1,1}.val-gammaln(alpha0) + gammaln(alpha0+N)            % final objective value
 
function nl = Gaussian_nll(X,mu_k,sigma_k,sigma_hat)        % Gaussian negative log-likelihood
   nl = (1/(2*(sigma_k+sigma_hat^2)))*sum((abs(X-mu_k)).^2)+(1/2)*log(sigma_k+sigma_hat^2); 
end

