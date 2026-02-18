function [mu,CJ,Tf_base,X0] = rs3_core_family_ic(name)
%RS3_CORE_FAMILY_IC  Baseline family initial conditions (reused core).
% Returns:
%   mu, CJ, Tf_base, X0 (reduced IC: [x;y;theta])
%
% Copied from baseline nested family_ic(name). New pipeline calls this as
% an allowed core reuse.
switch name
    case 'Lyapunov L1'
        mu=0.012150584270572; CJ=3.1293946655134501; Tf_base=2.9462815997965244;
        X04=[8.1152455023525571E-1; 0; 0; 2.5619581931804480E-1];
    case 'Lyapunov L2'
        mu=0.012150584270572; CJ=3.1300120370581501; Tf_base=3.4745552983343813;
        X04=[1.1010324052874807E+0; 0; 0; 2.6793038963934057E-1];
    case 'Cycler 21'
        mu=0.012150584270572; CJ=3.129389531054557; Tf_base=1.944007511640499e+01;
        X04=[7.237366530581342e-01; 0; 0; 4.137638707606867e-01];
    case 'Cycler 11a'
        mu=0.012150584270572; CJ=3.130000181904695; Tf_base=9.686582928190960;
        X04=[-0.810871347820831; 0; 0; -0.120835906310415];
    case 'Cycler 11b'
        mu=0.012150584270572; CJ=3.130001361366034; Tf_base=12.825036608914399;
        X04=[-0.759918239928426; 0; 0; -0.321726260888394];
    case 'Cycler 32'
        mu=0.012150584270572; CJ=3.130100479449716; Tf_base=18.044908875795237;
        X04=[-0.275704668101584; 0; 0; -2.112229772861937];
    case 'Resonant 2to1 Stable'
        mu=0.01215058427057155; CJ=3.1300000000000; Tf_base=6.0925;
        X04=[0.4493; 0; 0; 1.1826];
    case 'Resonant 2to1 Unstable'
        mu=0.01215058427057155; CJ=3.130000000000106; Tf_base=7.148969354042406;
        X04=[0.8695767910512225; 0; 0; -0.2689705299941745];
    case 'Resonant 3to1 Stable'
        mu=0.01215058427057155; CJ=3.1300000000000; Tf_base=6.2670;
        X04=[-0.8079; 0; 0; 0.1381];
    case 'Resonant 3to1 Unstable'
        mu=0.01215058427057155; CJ=3.129999999998302; Tf_base=6.454965220965995;
        X04=[0.8135643068820766; 0; 0; -0.2530482037736657];
    case 'Resonant 5to2 Stable'
        mu=0.01215058427057155; CJ=3.1300; Tf_base=12.6025;
        X04=[0.2283; 0; 0; 2.2737];
    case 'Resonant 5to2 Unstable'
        mu=0.01215058427057155; CJ=3.1300; Tf_base=12.9814554916569;
        X04=[-0.272350179729592; 0; 0; -2.13458916469565];
    case 'Distant Prograde Orbit'
        mu=0.01215058427057155; CJ=3.1300939875833800E+0; Tf_base=2.5551377141840610;
        X04=[1.0609987415775668E+0; 0; 0; 4.1093349889297970E-1];
    otherwise
        error('Unknown family "%s"', name);
end
theta0 = atan2(X04(4), X04(3));
X0 = [X04(1); X04(2); theta0];
end
