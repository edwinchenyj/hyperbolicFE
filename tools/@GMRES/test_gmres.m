function test_gmres()
    n = 120;    
    A = r(n);
    Av = @(v) A*v;
    
    obj = GMRES();
    obj.HMV = Av;
    e = [];
    for i = 1:100
        rng(10);
        b = rand(n,1);
        sol1 = obj.solve(b);
        sol2 = A\b;
        e = [e max(max(sol1 - sol2))];
    end
    plot(e)
end

function A = r(n)
rng(1);
A = rand(n);
end
function A = SPD(n)
% generating a SPD matrix;
    rng(1);
    A = rand(n);
    A = A + A';
    A = A + n*eye(n); % make sure it's diagonally dominating, for gershgorin circle theroem
end