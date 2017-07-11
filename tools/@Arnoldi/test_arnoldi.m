function test_arnoldi
    n = 120;    
    A = r(n);
    Av = @(v) A*v;
    
    obj.HMV = Av;
    e1 = [];
    e2 = [];
    for i = 1:100
        rng(10);
        b = rand(n,1);
        obj = Arnoldi(Av, b, 10);
        [Vmp1, HmBar, D] = obj.ConstructBasis;
        if D == length(b)+1
            Hm = HmBar(1:end-1,:);
            Vm = Vmp1(:,1:end-1);
        else
            Hm = HmBar(1:D,1:D);
            Vm = Vmp1(:,1:D);
        end
        Vm' * Vmp1
        e1 = [e1 norm(A * Vm - Vmp1 * HmBar)];
        e2 = [e2 norm(Vm' * A * Vm - Hm)];
    end
    hold on
    plot(e1)
    plot(e2)
end

function A = r(n)
rng(1);
A = rand(n);
end