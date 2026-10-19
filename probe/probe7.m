function probe7()
% Small, complete reference set used as the committed test fixtures.
outdir = '/home/hilal/dev/wsstpy/probe/data7';
if ~exist(outdir,'dir'); mkdir(outdir); end

N = 256; t = (0:N-1)';
rng(0);
S = struct();
S.tone    = cos(2*pi*0.05*t) + 0.5*cos(2*pi*0.20*t);
S.chirp   = cos(2*pi*(0.02*t + 0.4*t.^2/N));
S.noise   = randn(N,1);
S.short   = cos(2*pi*0.13*t) + 0.3*sin(2*pi*0.31*t);
S.impulse = double(t == floor(N/3));

% odd length
Nodd = 257; tod = (0:Nodd-1)';
S.odd = cos(2*pi*0.13*tod) + 0.3*sin(2*pi*0.31*tod);

W = struct();
for sig = {'tone','chirp','noise','short','impulse'}
    s = sig{1};
    [W.(['sst_' s]), W.(['f_' s])] = wsst(S.(s));
end
[W.sst_odd, W.f_odd] = wsst(S.odd);
[W.sst_vpo10, W.f_vpo10] = wsst(S.tone, 'VoicesPerOctave', 10);
[W.sst_extend, W.f_extend] = wsst(S.tone, 'ExtendSignal', true);
[W.sst_morse, W.f_morse]   = wsst(S.tone, 'Morse');
[W.sst_morse30, W.f_morse30] = wsst(S.tone, 'Morse', 'WaveletParameters', [3 30]);
[W.sst_bump, W.f_bump]     = wsst(S.tone, 'bump');
[W.sst_fs, W.f_fs]         = wsst(S.tone, 100);

C = struct();
C.x = S.tone;
[C.cfs_morse, C.f_morse] = cwt(C.x, 1);
C.icwt_full = icwt(C.cfs_morse);
C.icwt_band = icwt(C.cfs_morse, 'Morse', C.f_morse, [0.03 0.08]);
C.icwt_vpo  = icwt(C.cfs_morse, 'VoicesPerOctave', 10);
fb = cwtfilterbank('SignalLength', N, 'Wavelet', 'amor', ...
    'VoicesPerOctave', 32, 'SamplingFrequency', 1, 'Boundary', 'reflection');
C.cfs_amor = wt(fb, C.x);
C.f_amor = flip(fb.WaveletCenterFrequencies())';
C.icwt_amor = icwt(C.cfs_amor, 'amor');
[C.cfs_lim, C.f_lim] = cwt(C.x, 'FrequencyLimits', [0.02 0.25]);
C.icwt_lim = icwt(C.cfs_lim, 'Morse', C.f_lim, [0.03 0.1]);

I = struct();
I.full = iwsst(W.sst_tone);
I.band = iwsst(W.sst_tone, W.f_tone, [0.03 0.08]);
I.morse = iwsst(W.sst_morse, 'Morse');

G = struct();
for pen = [0 1 5 100]
    fld = sprintf('pen%d', round(pen));
    [G.(fld).fridge, G.(fld).iridge] = wsstridge(W.sst_tone, pen);
end
[G.r2.fridge, G.r2.iridge] = wsstridge(W.sst_tone, 5, 'NumRidges', 2);
[G.r3.fridge, G.r3.iridge] = wsstridge(W.sst_tone, 5, 'NumRidges', 3, 'NumFrequencyBins', 6);
[G.chirp.fridge, G.chirp.iridge] = wsstridge(W.sst_chirp, 2);
[G.noise.fridge, G.noise.iridge] = wsstridge(W.sst_noise, 3);
I.ridge = iwsst(W.sst_tone, G.r2.iridge);

save(fullfile(outdir,'ref7.mat'), 'S', 'W', 'C', 'I', 'G', '-v7');
fprintf('wsst %dx%d, cwt %dx%d, iridge %dx%d\n', size(W.sst_tone), size(C.cfs_morse), size(G.r2.iridge));
fprintf('DONE\n');
end
