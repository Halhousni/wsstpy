function probe2()
% Dump cwtfilterbank internals + wsst reference outputs for validation.
outdir = '/home/hilal/dev/wsstpy/probe/data';
if ~exist(outdir,'dir'); mkdir(outdir); end

cfgs = {256,10,'periodic'; 512,32,'periodic'; 1000,32,'periodic'; ...
        1024,32,'periodic'; 1024,32,'reflection'; 4096,48,'periodic'};
for k = 1:size(cfgs,1)
    N = cfgs{k,1}; vpo = cfgs{k,2}; bnd = cfgs{k,3};
    fb = cwtfilterbank('SignalLength',N,'VoicesPerOctave',vpo,'Wavelet','amor', ...
        'FrequencyLimits',[0 0.5],'Boundary',bnd,'SamplingFrequency',1);
    S = struct();
    S.N = N; S.vpo = vpo; S.bnd = bnd;
    S.scales = fb.Scales(:)';
    S.omega = fb.getFrequencyGrid();
    S.cf = fb.WaveletCenterFrequencies();
    S.psi = single(fb.PsiDFT);
    S.nyq = fb.NyquistBin;
    try, S.flim = fb.FrequencyLimits; catch, S.flim = [NaN NaN]; end
    save(fullfile(outdir,sprintf('fb_%d_%d_%s.mat',N,vpo,bnd)),'-struct','S','-v7');
    fprintf('N=%d vpo=%d %s: scales=%d nfft=%d flim=[%g %g] nyq=%d\n', N, vpo, bnd, ...
        numel(S.scales), numel(S.omega), S.flim(1), S.flim(2), S.nyq);
end

rng(0);
N = 1024;
sigs = {};
sigs{1} = cos(2*pi*0.1*(0:N-1)');
sigs{2} = cos(2*pi*0.05*(0:N-1)') + 0.5*cos(2*pi*0.2*(0:N-1)') + 0.1*randn(N,1);
sigs{3} = cos(2*pi*(0.02*(0:N-1)' + 0.4*(0:N-1)'.^2/N));
sigs{4} = randn(N,1);
sigs{5} = exp(-((0:N-1)'-N/2).^2/(2*(N/16)^2));
SSTs = cell(1,numel(sigs)); Fs = cell(1,numel(sigs));
SSTv = cell(1,numel(sigs)); Fv = cell(1,numel(sigs));
SSTr = cell(1,numel(sigs)); Fr = cell(1,numel(sigs));
for j = 1:numel(sigs)
    [SST,F] = wsst(sigs{j});                 SSTs{j} = single(SST); Fs{j} = F;
    [SST2,F2] = wsst(sigs{j},'VoicesPerOctave',10); SSTv{j} = single(SST2); Fv{j} = F2;
    [SST3,F3] = wsst(sigs{j},'ExtendSignal',true);  SSTr{j} = single(SST3); Fr{j} = F3;
end
save(fullfile(outdir,'wsst_ref.mat'),'sigs','SSTs','Fs','SSTv','Fv','SSTr','Fr','-v7');

for N2 = [37 100 1001]
    x = cos(2*pi*0.13*(0:N2-1)') + 0.3*sin(2*pi*0.31*(0:N2-1)');
    SST = wsst(x);
    [SST10,F10] = wsst(x,'VoicesPerOctave',10);
    eval(sprintf('SST_%d = single(SST); x_%d = x; F_%d = F10; SST10_%d = single(SST10);', N2, N2, N2, N2));
end
save(fullfile(outdir,'wsst_short.mat'),'SST_37','SST_100','SST_1001','x_37','x_100','x_1001', ...
    'F_37','F_100','F_1001','SST10_37','SST10_100','SST10_1001','-v7');
fprintf('DONE\n');
end
