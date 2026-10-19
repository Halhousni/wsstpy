function probe6()
% Per-signal references (chirp / noise / short) so every waveform is checked.
outdir = '/home/hilal/dev/wsstpy/probe/data6';
if ~exist(outdir,'dir'); mkdir(outdir); end

Ns = [256 512 1000 1024 1025 4096];
rng(0);
S = struct(); W = struct();
for k = 1:numel(Ns)
    N = Ns(k); t = (0:N-1)';
    S.(sprintf('chirp_%d', N)) = cos(2*pi*(0.02*t + 0.4*t.^2/N));
    S.(sprintf('noise_%d', N)) = randn(N,1);
    S.(sprintf('short_%d', N)) = cos(2*pi*0.13*t) + 0.3*sin(2*pi*0.31*t);
    S.(sprintf('impulse_%d', N)) = double(t == floor(N/3));
end

for k = 1:numel(Ns)
    N = Ns(k);
    for sig = {'chirp','noise','short','impulse'}
        s = sig{1};
        x = S.(sprintf('%s_%d', s, N));
        [sst, f] = wsst(x);
        W.(sprintf('sst_%s_%d', s, N)) = sst;
        W.(sprintf('f_%s_%d', s, N)) = f;
    end
end

% ridge on the chirp, from MATLAB's own SST
xc = S.chirp_1024;
sstc = W.sst_chirp_1024;
[G.chirp2.fridge, G.chirp2.iridge] = wsstridge(sstc, 2);
[G.chirp5.fridge, G.chirp5.iridge] = wsstridge(sstc, 5);
[G.chirp0.fridge, G.chirp0.iridge] = wsstridge(sstc, 0);
% noise case: many near-ties
sstn = W.sst_noise_1024;
[G.noise.fridge, G.noise.iridge] = wsstridge(sstn, 3);

save(fullfile(outdir,'ref6.mat'), 'S', 'W', 'G', '-v7');
fprintf('saved %d fields\n', numel(fieldnames(W)));
fprintf('DONE\n');
end
