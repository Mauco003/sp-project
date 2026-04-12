clear
close all
clc
%ciao a tutti
%% extract data

load("tracesbar1.mat");

varNames = ["pbar2438", "pbar2439", "pbar2440", "pbar2441", "pbar2442", ...
    "pbar2443", "pbar2444", "pbar2445", "pbar2446"];

p = cell(9);
for i = 1:9
    p{i} = eval(varNames(i));
    clear(varNames(i))
end

% correct mismatch in data
p7 = p{7};
temp = p7(:, 2);
p7(:, 2) = p7(:, 3);
p7(:, 3) = temp;

p{7} = p7;
clear temp p7

%% problem characteristic definition

Ts = 1/1e3; % sampling time

diamExt = 160;
diamInt = 100;
dx = (diamExt-diamInt)/2;
h = 290;

ksi = [0.68, 0.18, 0.14]; % mass percentage AP, Al, HTPB [-]
rho = [1.95, 2.7, 0.92];  % density [g/cm^3]
rho = rho/1e3;  % density [g/mm^3]

rhoP = 1/sum(ksi./rho);  % average density [g/mm^3]
VolTot = h*pi*(diamExt^2-diamInt^2)/4;  % total volume [mm^3]
MTot = VolTot * rhoP; % total mass [g]

diamThroat = [28.80, 25.26, 21.81]; % throat diameter [mm] (sorted by lower to higher pressure)
Athroat = pi*diamThroat.^2/4; % troat area [mm^2]


%%
rb = zeros(3, 9);
pEff = zeros(3, 9);
cStar = zeros(3, 9);

for i = 1:9
    pI = p{i};

    t = Ts*(1:size(pI, 1));
    pLow = pI(:, 1);
    pMid = pI(:, 2);
    pHigh = pI(:, 3);

    [tALow, tDLow, PALow, PDLow] = actionTime(pLow, Ts, 5);
    [tAMid, tDMid, PAMid, PDMid] = actionTime(pMid, Ts, 5);
    [tAHigh, tDHigh, PAHigh, PDHigh] = actionTime(pHigh, Ts, 5);

    actionTimeIdx = t>=tALow & t <=tDHigh;
    pRefLow = trapz(t(t>=tALow & t <=tDLow), pLow(t>=tALow & t <=tDLow))/(2*(tDLow - tALow));
    pRefMid = trapz(t(t>=tAMid & t <=tDMid), pMid(t>=tAMid & t <=tDMid))/(2*(tDMid - tAMid));
    pRefHigh = trapz(t(t>=tAHigh & t <=tDHigh), pHigh(t>=tAHigh & t <=tDHigh))/(2*(tDHigh - tAHigh));

    [tBLow, tCLow, percBurnLow] = burningTime(pLow, Ts, pRefLow);
    [tBMid, tCMid, percBurnMid] = burningTime(pMid, Ts, pRefMid);
    [tBHigh, tCHigh, percBurnHigh] = burningTime(pHigh, Ts, pRefHigh);

    rb(1, i) = dx/(tCLow-tBLow);
    rb(2, i) = dx/(tCMid-tBMid);
    rb(3, i) = dx/(tCHigh-tBHigh);

    pEff(1, i) = trapz(t(t>=tBLow & t <=tCLow), pLow(t>=tBLow & t <=tCLow, 1))/(tCLow - tBLow);
    pEff(2, i) = trapz(t(t>=tBMid & t <=tCMid), pMid(t>=tBMid & t <=tCMid, 1))/(tCMid - tBMid);
    pEff(3, i) = trapz(t(t>=tBHigh & t <=tCHigh), pHigh(t>=tBHigh & t <=tCHigh, 1))/(tCHigh - tBHigh);

    cStar(1, i) = trapz(t(t>=tBLow & t <=tCLow), pLow(t>=tBLow & t <=tCLow, 1)*Athroat(1))/MTot;
    cStar(2, i) = trapz(t(t>=tBMid & t <=tCMid), pMid(t>=tBMid & t <=tCMid, 1)*Athroat(2))/MTot;
    cStar(3, i) = trapz(t(t>=tBHigh & t <=tCHigh), pHigh(t>=tBHigh & t <=tCHigh, 1)*Athroat(3))/MTot;
end


[a, aSigma, n, nSigma, R2] = Uncertainty(pEff(:), rb(:));

cStar = cStar * 1e2; % convert from [bar mm^2 s / g] to [m / s]
[cStarSigma, cStarAvg] = std(cStar(:));


%%
close all
[tPredictLow, pPredictLow, rbPredictLow] = balisticPredict(a, n, rhoP, cStarAvg, diamExt, diamInt, h, Athroat(1));
[tPredictMid, pPredictMid, rbPredictMid] = balisticPredict(a, n, rhoP, cStarAvg, diamExt, diamInt, h, Athroat(2));
[tPredictHigh, pPredictHigh, rbPredictHigh] = balisticPredict(a, n, rhoP, cStarAvg, diamExt, diamInt, h, Athroat(3));

%%
close all
for i = 1:9
    pI = p{i};
    t = Ts*(1:size(pI, 1));
    for j = 1:3
        figure(j)
        hold on
        plot(t, pI(:, j))
        hold off
    end
end

figure(1)
hold on
plot(tPredictLow, pPredictLow, "r", "LineWidth", 3)
hold off

figure(2)
hold on
plot(tPredictMid, pPredictMid, "r", "LineWidth", 3)
hold off


figure(3)
hold on
plot(tPredictHigh, pPredictHigh, "r", "LineWidth", 3)
hold off

for i = 1:3
    figure(i)
    legend([varNames, "predicted"], "FontSize", 12, "FontName", 'Times New Roman')
    ax = gca();
    ax.FontSize = 12;
    ylabel("P [bar]", "FontSize", 12, "FontName",'Times New Roman')
    xlabel("t [s]", "FontSize", 12, "FontName", 'Times New Roman')
    xlim([0, 5])
    box on
    grid on
end


%%

if false
    hFig = figure;
    set(hFig, "windowState", "maximized")
    
    for i = 1:9
        
        subplot(3, 3, i);
        pI = p{i};
        
        t = Ts*(1:size(pI, 1));
        pLow = pI(:, 1);
        pMid = pI(:, 2);
        pHigh = pI(:, 3);
    
        hLow = plot(t, pLow);
        hold on
        hMid = plot(t, pMid);
        hHigh = plot(t, pHigh);
        
        xlim([0, 5.2])
        ylim([0, 90])
    
        [tALow, tDLow, PALow, PDLow] = actionTime(pLow, Ts, 5);
        [tAMid, tDMid, PAMid, PDMid] = actionTime(pMid, Ts, 5);
        [tAHigh, tDHigh, PAHigh, PDHigh] = actionTime(pHigh, Ts, 5);
        
        xline([tALow, tDLow], "--", "Color", hLow.Color)
        xline([tAMid, tDMid], "--", "Color", hMid.Color)
        xline([tAHigh, tDHigh], "--", "Color", hHigh.Color)
    
        yline([PALow, PDLow], "--", "Color", hLow.Color)
        yline([PAMid, PDMid], "--", "Color", hMid.Color)
        yline([PAHigh, PDHigh], "--", "Color", hHigh.Color)
    
        legend([hLow, hMid, hHigh], "low", "medium", "high")
    end
end

%%


function [maxP, idx] = findMax(p)
    pMod = p;
    pMod(1:200) = 0;
    [maxP, idx] = max(pMod);
end


function [tA, tD, PA, PD] = actionTime(p, Ts, perc)
    [maxP, idxPMax] = findMax(p);

    idxPA = find(p(idxPMax:-1:1)/maxP < perc/100, 1, "first");
    idxPA = idxPMax - idxPA + 1;

    idxPD = find(p(idxPMax:end)/maxP < perc/100, 1, "first");
    idxPD = idxPMax + idxPD - 1;

    PA = p(idxPA);
    PD = p(idxPD);
    tA = idxPA*Ts;
    tD = idxPD*Ts;
end

function [tB, tC, perc] = burningTime(p, Ts, pRef)
    [maxP, idxPMax] = findMax(p);

    perc = pRef/maxP*100;

    idxPB = find(p(idxPMax:-1:1) < pRef, 1, "first");
    idxPB = idxPMax - idxPB + 1;

    idxPC = find(p(idxPMax:end) < pRef, 1, "first");
    idxPC = idxPMax + idxPC - 1;

    tB = idxPB*Ts;
    tC = idxPC*Ts;
    
end


function [t, p, rb] = balisticPredict(a, n, rhoP, cStar, diamExt, diamInt, h, Athroat)
    
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