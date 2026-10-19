function probe5()
% Double-precision references for the validation report.
outdir = '/home/hilal/dev/wsstpy/probe/data5';
if ~exist(outdir,'dir'); mkdir(outdir); end

R = struct();
R.matlab_version = version;
R.date = datestr(now, 31);

% --- signals -------------------------------------------------------------
Ns = [256 512 1000 1024 1025 4096];
rng(0);
S = struct();
for k = 1:numel(Ns)
    N = Ns(k); t = (0:N-1)';
    S.(sprintf('tone_%d', N))    = cos(2*pi*0.05*t) + 0.5*cos(2*pi*0.20*t);
    S.(sprintf('chirp_%d', N))   = cos(2*pi*(0.02*t + 0.4*t.^2/N));
    S.(sprintf('noise_%d', N))   = randn(N,1);
    S.(sprintf('short_%d', N))   = cos(2*pi*0.13*t) + 0.3*sin(2*pi*0.31*t);
end

% --- wsst in double, all wavelets ---------------------------------------
W = struct();
for k = 1:numel(Ns)
    N = Ns(k);
    fld = sprintf('tone_%d', N);
    x = S.(fld);
    [W.(sprintf('sst_%d', N)), W.(sprintf('f_%d', N))] = wsst(x);
    [W.(sprintf('sst_vpo10_%d', N)), W.(sprintf('f_vpo10_%d', N))] = ...
        wsst(x, 'VoicesPerOctave', 10);
    [W.(sprintf('sst_ext_%d', N)), W.(sprintf('f_ext_%d', N))] = ...
        wsst(x, 'ExtendSignal', true);
end
x = S.tone_1024;
[W.sst_morse, W.f_morse]   = wsst(x, 'Morse');
[W.sst_morse30, W.f_morse30] = wsst(x, 'Morse', 'WaveletParameters', [3 30]);
[W.sst_morse_g2, W.f_morse_g2] = wsst(x, 'Morse', 'WaveletParameters', [2 40]);
[W.sst_bump, W.f_bump]     = wsst(x, 'bump');
[W.sst_fs, W.f_fs]         = wsst(x, 100);

% --- cwt / icwt ---------------------------------------------------------
C = struct();
[C.cfs_morse, C.f_morse] = cwt(x, 1);
C.icwt_full = icwt(C.cfs_morse);
C.icwt_band = icwt(C.cfs_morse, 'Morse', C.f_morse, [0.03 0.08]);
C.icwt_vpo  = icwt(C.cfs_morse, 'VoicesPerOctave', 10);
fb = cwtfilterbank('SignalLength', 1024, 'Wavelet', 'amor', ...
    'VoicesPerOctave', 32, 'SamplingFrequency', 1, 'Boundary', 'reflection');
C.cfs_amor = wt(fb, x);
C.f_amor = flip(fb.WaveletCenterFrequencies())';
C.icwt_amor = icwt(C.cfs_amor, 'amor');
[C.cfs_lim, C.f_lim] = cwt(x, 'FrequencyLimits', [0.02 0.25]);
C.icwt_lim = icwt(C.cfs_lim, 'Morse', C.f_lim, [0.03 0.1]);

% --- iwsst --------------------------------------------------------------
I = struct();
sst = W.sst_1024;
I.full = iwsst(sst);
I.band = iwsst(sst, W.f_1024, [0.03 0.08]);
[~, ir2] = wsstridge(sst, 5, 'NumRidges', 2);
I.ridge = iwsst(sst, ir2);
I.morse = iwsst(W.sst_morse, 'Morse');
I.amor2 = iwsst(sst, 'amor');

% --- ridges -------------------------------------------------------------
G = struct();
x3 = S.tone_1024;
sst3 = W.sst_1024;
for pen = [0 0.5 1 5 10 100]
    fld = sprintf('pen%d', round(pen*10));
    [G.(fld).fridge, G.(fld).iridge] = wsstridge(sst3, pen);
end
[G.r2.fridge, G.r2.iridge] = wsstridge(sst3, 5, 'NumRidges', 2);
[G.r3.fridge, G.r3.iridge] = wsstridge(sst3, 5, 'NumRidges', 3, 'NumFrequencyBins', 6);
xc = S.chirp_1024;
[G.chirp.fridge, G.chirp.iridge] = wsstridge(wsst(xc), 2);

save(fullfile(outdir, 'ref5.mat'), 'R', 'S', 'W', 'C', 'I', 'G', '-v7');
fprintf('saved: wsst %dx%d, cwt %dx%d, icwt %d\n', ...
    size(W.sst_1024), size(C.cfs_morse), numel(C.icwt_full));
fprintf('DONE\n');
end
