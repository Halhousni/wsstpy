running wsst comparisons ...
running cwt / icwt comparisons ...
running iwsst comparisons ...
running ridge comparisons ...

MATLAB 25.2.0.3177638 (R2025b) Update 5 | Python 3.13.5 | numpy 2.5.4 | scipy 1.18.1

| comparison | max abs diff | max rel diff | RMS diff |
| --- | --- | --- | --- |
| wsst tone N=256 | 2.168e-15 | 1.008e-14 | 7.713e-17 |
| wsst vpo=10 N=256 | 2.959e-16 | 1.372e-15 | 2.455e-17 |
| wsst extend N=256 | 2.445e-15 | 1.137e-14 | 7.268e-17 |
| wsst tone N=512 | 2.478e-15 | 1.152e-14 | 6.406e-17 |
| wsst vpo=10 N=512 | 4.282e-16 | 1.981e-15 | 2.507e-17 |
| wsst extend N=512 | 2.162e-15 | 1.005e-14 | 4.764e-17 |
| wsst tone N=1000 | 7.005e-16 | 3.256e-15 | 3.716e-17 |
| wsst vpo=10 N=1000 | 2.243e-16 | 1.042e-15 | 1.438e-17 |
| wsst extend N=1000 | 1.670e-15 | 7.761e-15 | 3.907e-17 |
| wsst tone N=1024 | 2.205e-15 | 1.025e-14 | 4.598e-17 |
| wsst vpo=10 N=1024 | 4.475e-16 | 2.064e-15 | 1.663e-17 |
| wsst extend N=1024 | 2.071e-15 | 9.628e-15 | 4.159e-17 |
| wsst tone N=1025 | 2.134e-15 | 9.919e-15 | 4.911e-17 |
| wsst vpo=10 N=1025 | 4.679e-16 | 2.160e-15 | 2.481e-17 |
| wsst extend N=1025 | 2.103e-15 | 9.776e-15 | 4.600e-17 |
| wsst tone N=4096 | 2.168e-15 | 1.008e-14 | 3.399e-17 |
| wsst vpo=10 N=4096 | 3.310e-16 | 1.535e-15 | 1.402e-17 |
| wsst extend N=4096 | 2.448e-15 | 1.138e-14 | 3.464e-17 |
| wsst morse (3,60) | 2.178e-15 | 1.330e-14 | 5.159e-17 |
| wsst morse (3,30) | 2.185e-15 | 9.321e-15 | 5.098e-17 |
| wsst morse (2,40) | 2.567e-15 | 1.285e-14 | 8.927e-17 |
| wsst bump | 1.939e-15 | 2.662e-14 | 7.721e-17 |
| wsst fs=100 | 2.205e-15 | 1.025e-14 | 4.598e-17 |
| cwt Morse default | 8.111e-15 | 8.136e-15 | 1.544e-15 |
| cwt Morse limits [0.02 0.25] | 5.997e-15 | 6.035e-15 | 1.341e-15 |
| cwt amor 32 voices | 2.143e-14 | 2.085e-14 | 4.735e-15 |
| icwt full | 6.217e-15 | 4.129e-15 | 2.823e-15 |
| icwt band [0.03 0.08] | 4.663e-15 | 4.313e-15 | 2.599e-15 |
| icwt vpo=10 | 6.217e-15 | 4.129e-15 | 2.823e-15 |
| icwt amor | 3.461e-07 | 7.193e-08 | 1.819e-07 |
| icwt Morse limits [0.03 0.1] | 4.663e-15 | 4.279e-15 | 2.525e-15 |
| iwsst full | 3.651e-08 | 2.429e-08 | 1.922e-08 |
| iwsst band [0.03 0.08] | 2.583e-08 | 2.429e-08 | 1.716e-08 |
| iwsst by ridge | 2.545e-08 | 2.429e-08 | 1.356e-08 |
| iwsst Morse | 4.441e-15 | 2.952e-15 | 2.083e-15 |
| wsstridge penalty=0 | 1024/1024 identical | 0 |  |
| wsstridge penalty=0.5 | 1024/1024 identical | 0 |  |
| wsstridge penalty=1 | 1024/1024 identical | 0 |  |
| wsstridge penalty=5 | 1024/1024 identical | 0 |  |
| wsstridge penalty=10 | 1024/1024 identical | 0 |  |
| wsstridge penalty=100 | 1024/1024 identical | 0 |  |
| wsstridge 2 ridges | 2048/2048 identical | 0 |  |
| wsstridge 3 ridges, 6 bins | 3072/3072 identical | 0 |  |
| wsstridge chirp penalty=2 | 1024/1024 identical | 0 |  |

rows: 63
