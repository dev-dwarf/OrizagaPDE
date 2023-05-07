run('Configuration.m');
generateData = true;

for equation=HEAT:ALLENCAHN
    for dxMagnitude=0:dxMagnitudes
        dtMagnitude = 1;
        N = 25*2^dxMagnitude;
        run("NoFlux.m");
        run("Periodic.m");
    end

    for dtMagnitude = 0:dtMagnitudes
        dxMagnitude = 1;
        N = 25*2^dxMagnitude;
        run("NoFlux.m");
        run("Periodic.m");
    end
end

clear generateData;

