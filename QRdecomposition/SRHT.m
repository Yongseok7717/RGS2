%% Subsampled randomized Hadamard Transform (SRHT).
function Theta = SRHT(n,k)
    D = randi([0 1], n,1)*2 - 1;
    N = 2^ceil(log(n)/log(2));
    perm = randperm(N,k);
    select = @(t,ind) t(ind);
    Theta = @(t) (1/sqrt(k)) * select(myfwht(D.*t),perm);
end


%% Fast Walsh Hadamard Transform
function z = myfwht(a)
    h = 1;
    n = length(a);
    N = 2^ceil(log(n)/log(2));
    z = zeros(N,1);
    z(1:n) = a;
    while h < N
        for i = 1:2*h:N
            for j = i:(i + h - 1)
                x = z(j);
                y = z(j + h);
                z(j) = x + y;
                z(j + h) = x - y;
            end
        end
        h = 2*h;
    end
end