function [t, p, rb] = computeBurn(a, n, rhoP, cStar, diamExt, diamInt, h, Athroat)
    
    rExt = diamExt/2;
    rInt0 = diamInt/2;

    x0 = [rInt0; h];
    cStar = cStar/1e2;

    SRM.a = a;
    SRM.rhoP = rhoP;
    SRM.cStar = cStar;
    SRM.rExt = rExt;
    SRM.Athroat = Athroat;
    SRM.n = n;

    options = odeset("Events", @(t, x) eventFunc(t, x, rExt), "RelTol", 1e-9, "AbsTol", 1e-10);

    [t, x] = ode45(@(t, x) burnODE(t, x, SRM), [0, inf], x0, options);
    rInt = x(:, 1)';
    h = x(:, 2)';
    Ab = 2*pi*(rExt^2-rInt.^2)+2*pi*h.*rInt;
    p = (a.*rhoP.*cStar.*Ab./Athroat).^(1./(1-n));
    rb = a.*(p.^n);
end

function dx = burnODE(~, x, SRM)
    a = SRM.a;
    rhoP = SRM.rhoP;
    cStar = SRM.cStar; % recovering unit measure balanced to other terms
    rExt = SRM.rExt;
    Athroat = SRM.Athroat;
    n = SRM.n;
    
    rInt = x(1);
    h = x(2);
    Ab = 2*pi*(rExt^2-rInt.^2)+ 2*pi*h.*rInt;
    p = (a.*rhoP.*cStar.*Ab./Athroat).^(1./(1-n));
    rb = a.*(p.^n);

    dx = zeros(2, 1);
    dx(1) = +rb;
    dx(2) = -2*rb;
end

function [value,isterminal,direction] = eventFunc(~,x, rExt)
    value = ~(rExt - x(1)<1e-14 || x(2)<1e-14);
    isterminal = 1;
    direction = 0;
end