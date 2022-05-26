clearvars
clc

AT = ActinTracker;

process(AT, ...
    'D:\Projects\ALMC Tickets\T50-Lehman-Actin\data\T370 Fu Data\BetaWT_50_1_002.nd2', ...
    'D:\Projects\ALMC Tickets\T50-Lehman-Actin\results')

%%

AD = ActinData;

AD = importdata(AD, 'D:\Projects\ALMC Tickets\T50-Lehman-Actin\results\BetaWT_50_1_002.mat');
AD = analyze(AD);

imshow(showlabels(AD, 5))
