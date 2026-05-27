function [mu,CJ,Tf_base,X0] = cr3bp_family_ic(name)
%CR3BP_FAMILY_IC  Initial conditions for the 13 representative periodic orbits.
%
% These are the orbits listed in Table 2 of the paper, all corrected to the
% common Jacobi constant CJ = 3.1294 in the planar Earth-Moon CR3BP.
% Initial conditions were computed with correct_po_to_cj313_v2.m and
% verified against the instability rates and periods in Table 2.
%
% Usage:
%   [mu, CJ, Tf_base, X0] = cr3bp_family_ic(name)
%
% Input:
%   name    — family name string (see cases below)
%
% Outputs:
%   mu      — Earth-Moon mass parameter (0.012150584270572)
%   CJ      — Jacobi constant (nominally 3.1294; see note below)
%   Tf_base — nondimensional period (one full orbit, TU units)
%   X0      — reduced 3D IC: [x; y; theta] on the y=0 section
%
% Family name → paper abbreviation (Table 2) mapping:
%   'Lyapunov L1'           → LL1   (12.811 days)
%   'Lyapunov L2'           → LL2   (15.117 days)
%   'Cycler 11a'            → C11a  (42.140 days)
%   'Cycler 11b'            → C11b  (55.995 days)
%   'Cycler 21'             → C21   (84.533 days)
%   'Cycler 32'             → C32   (78.613 days)
%   'Resonant 2to1 Stable'  → R21-S (26.500 days)
%   'Resonant 2to1 Unstable'→ R21-U (31.039 days)
%   'Resonant 3to1 Stable'  → R31-S (27.252 days)
%   'Resonant 3to1 Unstable'→ R31-U (28.066 days)
%   'Resonant 5to2 Stable'  → R52-S (54.802 days)
%   'Resonant 5to2 Unstable'→ R52-U (56.436 days)
%   'Distant Prograde Orbit'→ DPO   (11.184 days)
%
% Note: a handful of families could not be converged to CJ = 3.1294 exactly;
% their stored CJ values differ in the last few digits (e.g., Cycler 21 at
% CJ = 3.129389...). This is noted in Table 2 of the paper and does not
% affect the shared-energy-manifold analysis because all representatives
% lie within < 0.001% of the nominal value.
switch name
    case 'Lyapunov L1'
        mu=0.012150584270572; CJ=3.1294; Tf_base=2.946253150022597;
        X04=[0.8115256290557147; 0; 0; 0.2561843220006502];
    case 'Lyapunov L2'
        mu=0.012150584270572; CJ=3.1294; Tf_base=3.476412471522109;
        X04=[1.100554841329441; 0; 0; 0.2702199577829618];
    case 'Cycler 21'
        mu=0.012150584270572; CJ=3.129389531054557; Tf_base=1.944007511640499e+01;
        X04=[7.237366530581342e-01; 0; 0; 4.137638707606867e-01];
    case 'Cycler 11a'
        mu=0.012150584270572; CJ=3.1294; Tf_base=9.691077443795967;
        X04=[-0.8116406668238326; 0; 0; -0.1185905575975559];
    case 'Cycler 11b'
        mu=0.012150584270572; CJ=3.1294; Tf_base=12.87711489493616;
        X04=[-0.7601803401119139; 0; 0; -0.3218380227492717];
    case 'Cycler 32'
        mu=0.012150584270572; CJ=3.129399999999999; Tf_base=18.07860804203491;
        X04=[-0.2752114538522783; 0; 0; -2.115657313283594];
    case 'Resonant 2to1 Stable'
        mu=0.0121505842705716; CJ=3.1294; Tf_base=6.09412703751587;
        X04=[0.4485378415692721; 0; 0; 1.185506250841159];
    case 'Resonant 2to1 Unstable'
        mu=0.0121505842705716; CJ=3.1294; Tf_base=7.138022234549252;
        X04=[0.86975081578361; 0; 0; -0.2703857273991174];
    case 'Resonant 3to1 Stable'
        mu=0.0121505842705716; CJ=3.1294; Tf_base=6.267093371035164;
        X04=[-0.8081272737956029; 0; 0; 0.1389495551262821];
    case 'Resonant 3to1 Unstable'
        mu=0.0121505842705716; CJ=3.1294; Tf_base=6.454444248544365;
        X04=[0.8137553067908899; 0; 0; -0.2540548036743875];
    case 'Resonant 5to2 Stable'
        mu=0.0121505842705716; CJ=3.1294; Tf_base=12.60275494862142;
        X04=[0.2278881086717652; 0; 0; 2.277116891512287];
    case 'Resonant 5to2 Unstable'
        mu=0.0121505842705716; CJ=3.1294; Tf_base=12.97870170292832;
        X04=[-0.2719428329684943; 0; 0; -2.137466033981196];
    case 'Distant Prograde Orbit'
        mu=0.0121505842705716; CJ=3.1294; Tf_base=2.572014144294218;
        X04=[1.060820806800189; 0; 0; 0.4126719472859519];
    otherwise
        error('Unknown family "%s"', name);
end
theta0 = atan2(X04(4), X04(3));
X0 = [X04(1); X04(2); theta0];
end