function probe4()
outdir = '/home/hilal/dev/wsstpy/probe/data';
N = 1024; t = (0:N-1)';
x1 = cos(2*pi*0.05*t) + cos(2*pi*0.20*t);

% icwt with defaults (Morse, vpo 10)
[cfs, fc] = cwt(x1, 1);
xi_full = icwt(cfs);
xi_band = icwt(cfs, 'Morse', fc, [0.03 0.08]);
xi_vpo = icwt(cfs, 'VoicesPerOctave', 10);

% icwt for the Morlet filter bank
fb = cwtfilterbank('SignalLength',N,'Wavelet','amor','VoicesPerOctave',32, ...
    'SamplingFrequency',1,'Boundary','reflection');
cfs_amor = wt(fb, x1);
xi_amor = icwt(cfs_amor, 'Wavelet', 'amor');

% frequencies in Hz
[SSTfs, Ffs] = wsst(x1, 100);
% iwsst with Morse
[SSTm, Fm] = wsst(x1, 'Morse');
xrm = iwsst(SSTm, 'Morse');

save(fullfile(outdir,'wsst_ref3.mat'), 'x1','cfs','fc','xi_full','xi_band','xi_vpo', ...
    'cfs_amor','xi_amor','SSTfs','Ffs','SSTm','Fm','xrm','-v7');
fprintf('xi_full %g xi_amor %g\n', max(abs(xi_full)), max(abs(xi_amor)));
fprintf('DONE\n');
end
