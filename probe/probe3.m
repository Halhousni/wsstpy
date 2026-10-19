function probe3()
outdir = '/home/hilal/dev/wsstpy/probe/data';
N = 1024; t = (0:N-1)';
x1 = cos(2*pi*0.05*t) + cos(2*pi*0.20*t);
x2 = cos(2*pi*(0.05*t + 0.15*t.^2/N));      % chirp
x3 = cos(2*pi*0.05*t) + 0.7*cos(2*pi*0.20*t) + 0.3*cos(2*pi*0.38*t);
[SST1,F1] = wsst(x1);
[SST2,F2] = wsst(x2);
[SST3,F3] = wsst(x3);

% ridge extraction, several penalties / ridge counts
[fr_p0, ir_p0]   = wsstridge(SST1, 0);
[fr_p1, ir_p1]   = wsstridge(SST1, 1);
[fr_p10, ir_p10] = wsstridge(SST1, 10);
[fr_p100, ir_p100] = wsstridge(SST1, 100);
[fr_p5_2, ir_p5_2] = wsstridge(SST1, 5, 'NumRidges', 2);
[fr_p5_3b, ir_p5_3b] = wsstridge(SST3, 5, 'NumRidges', 3, 'NumFrequencyBins', 4);
[fr_ch, ir_ch] = wsstridge(SST2, 2);

% inverse
xr_full = iwsst(SST1);
xr_band = iwsst(SST1, F1, [0.03 0.08]);
xr_ridge = iwsst(SST1, ir_p5_2);

% cwt
cfs = cwt(x1, 1);
fb = cwtfilterbank('SignalLength',N,'Wavelet','amor','VoicesPerOctave',32,'SamplingFrequency',1);
cfs_def = wt(fb,x1);
[cfm, fcm] = cwt(x1);

% morse and bump wsst
[SSTm, Fm] = wsst(x1, 'Morse');
[SSTb, Fb] = wsst(x1, 'bump');
[SSTm2, Fm2] = wsst(x1, 'Morse', 'WaveletParameters', [3 30]);
[SSTm10, Fm10] = wsst(x1, 'Morse', 'VoicesPerOctave', 10);

save(fullfile(outdir,'wsst_ref2.mat'), 'x1','x2','x3','SST1','F1','SST2','F2','SST3','F3', ...
    'fr_p0','ir_p0','fr_p1','ir_p1','fr_p10','ir_p10','fr_p100','ir_p100', ...
    'fr_p5_2','ir_p5_2','fr_p5_3b','ir_p5_3b','fr_ch','ir_ch', ...
    'xr_full','xr_band','xr_ridge', ...
    'cfs','cfs_def','cfm','fcm', ...
    'SSTm','Fm','SSTb','Fb','SSTm2','Fm2','SSTm10','Fm10', '-v7');
fprintf('SST1 %dx%d, F1 %d, ir_p5_2 %dx%d\n', size(SST1), numel(F1), size(ir_p5_2));
fprintf('DONE\n');
end
