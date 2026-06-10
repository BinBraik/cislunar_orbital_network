function short_names = net_family_short_names()
%NET_FAMILY_SHORT_NAMES  Short label for each of the 13 periodic-orbit families.
%
%   Indices match the exact storage order in the sweep .mat file.
%
%   Returns
%     short_names : {13×1 cell} of character vectors

short_names = { ...
    'LyapL1';   ...  %  1  Lyapunov L1
    'LyapL2';   ...  %  2  Lyapunov L2
    'Cyc21';    ...  %  3  Cycler 21
    'Cyc11a';   ...  %  4  Cycler 11a
    'Cyc11b';   ...  %  5  Cycler 11b
    'Cyc32';    ...  %  6  Cycler 32
    'R21S';     ...  %  7  Resonant 2to1 Stable
    'R21U';     ...  %  8  Resonant 2to1 Unstable
    'R31S';     ...  %  9  Resonant 3to1 Stable
    'R31U';     ...  % 10  Resonant 3to1 Unstable
    'R52S';     ...  % 11  Resonant 5to2 Stable
    'R52U';     ...  % 12  Resonant 5to2 Unstable
    'DPO';      ...  % 13  Distant Prograde Orbit
    };

end
