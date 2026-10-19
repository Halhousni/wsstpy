% probe1: dump cwtfilterbank internals for wsst's exact configuration
N = 1024; vpo = 32;
fb = cwtfilterbank('SignalLength',N,'VoicesPerOctave',vpo,'Wavelet','amor', ...
    'FrequencyLimits',[0 0.5],'Boundary','periodic','SamplingFrequency',1);
fprintf('num scales: %d\n', numel(fb.Scales));
fprintf('scales(1:3): %.17g %.17g %.17g\n', fb.Scales(1),fb.Scales(2),fb.Scales(3));
fprintf('scales(end-2:end): %.17g %.17g %.17g\n', fb.Scales(end-2),fb.Scales(end-1),fb.Scales(end));
fprintf('PsiDFT size: %d x %d, class %s\n', size(fb.PsiDFT,1), size(fb.PsiDFT,2), class(fb.PsiDFT));
fprintf('sigpad: %d  NyquistBin: %d\n', fb.SignalPad, fb.NyquistBin);
cf = fb.WaveletCenterFrequencies();
fprintf('cf(1:3): %.17g %.17g %.17g  cf(end): %.17g\n', cf(1),cf(2),cf(3),cf(end));
om = fb.getFrequencyGrid();
fprintf('omega(1:4): %.17g %.17g %.17g %.17g  numel %d\n', om(1),om(2),om(3),om(4),numel(om));
% peak of row 1 and row end
[~,i1] = max(abs(fb.PsiDFT(1,:))); [~,i2] = max(abs(fb.PsiDFT(end,:)));
fprintf('row1 peak at omega=%.17g val=%.17g ; rowEnd peak at omega=%.17g val=%.17g\n', om(i1), fb.PsiDFT(1,i1), om(i2), fb.PsiDFT(end,i2));
fprintf('row1 nonzero count %d, rowEnd nonzero count %d\n', nnz(fb.PsiDFT(1,:)), nnz(fb.PsiDFT(end,:)));
fprintf('isreal PsiDFT: %d, any complex: %d\n', isreal(fb.PsiDFT), ~isreal(fb.PsiDFT));
save('/home/hilal/dev/wsstpy/probe/fb1024.mat','-v7.3');
