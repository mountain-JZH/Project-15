close all
clear all


HdN01 = filterN01;

FS_baseband = 3.2e9; % 3.2GHz
fs_baseband = 64e6; % 64MHz

%% [OFDM Parameter]

M = 16; % Modulation order for 16QAM 
nfft = 64; % IFFT
cplen = 0;
nSym = 1000; % Number of symbols 
numbSubcarriers = 16; % Number of data carrying subcarriers 
scs = 1e6; % 1MHz

nullIdx = [1:nfft/2 - numbSubcarriers/2, nfft/2 + numbSubcarriers/2 + 1 : nfft]';
numDataCarrs = nfft - length(nullIdx);

inSig = randi([0 M-1], numDataCarrs, nSym);
qamSym = qammod(inSig, M, 'UnitAveragePower', false);

% [Plot] [1] one OFDM symbol in F domain
firstSym = qamSym(:, 1);
X = zeros(nfft, 1);
dataIdx = setdiff(1:nfft, nullIdx);
X(dataIdx) = firstSym;

figure(1)
stem(0:nfft-1, abs(X), 'filled');
xlabel('Subcarrier Index');
ylabel('Magnitude');
title(' one OFDM symbol in F domain');
grid on;


outSig = ofdmmod(qamSym, nfft, cplen, nullIdx); % OFDM signal in T domain

%% [Limiter/Clipping]

P_avg = mean(abs(outSig).^2);      
RMS = sqrt(P_avg);                 
clip_ratio = 2.0;                  
limit = clip_ratio * RMS;          

outSig_clipped = outSig; % OFDM signal in T domain after clipping        

for k = 1:length(outSig)
    if abs(outSig(k)) > limit
        outSig_clipped(k) = outSig(k) / abs(outSig(k)) * limit;
    end
end

% [Plot] [2] OFDM before/after Clipping in T domain
figure(2); 

% [Plot] [2.1] before Clipping
subplot(1, 2, 1);
plot(abs(outSig));
ylim([0, 1]);
title('OFDM before Clipping in T domain');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% [Plot] [2.2] after Clipping  
subplot(1, 2, 2); 
plot(abs(outSig_clipped));
ylim([0, 1]);
title('OFDM after Clipping in T domain');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% [Plot] [3] PSD before & after Clipping
figure(3);

% [Plot] [3.1] PSD before Clipping
subplot(1, 3, 1);
pwelch(outSig, [], [], [], fs_baseband, 'centered');
title('PSD of OFDM Signal (Before Clipping)');
ylim([-150, -80]); 
grid on;

% [Plot] [3.2] PSD after Clipping
subplot(1, 3, 2);
pwelch(outSig_clipped, [], [], [], fs_baseband, 'centered');
title('PSD of OFDM Signal (After Clipping)');
ylim([-150, -80]); 
grid on;

%% [(c) = (a) - (b) full distortion]

% [Plot] [4] full distortion
distortion = outSig - outSig_clipped;                 % c = a - b
out_of_band_distortion = filter(HdN01, distortion);   % d


gdN01 = round(mean(grpdelay(HdN01,1024))); % 59


a_aligned = [zeros(gdN01,1); outSig];
a_aligned = a_aligned(1:length(out_of_band_distortion));
outSig_PAPR_reduction = a_aligned - out_of_band_distortion;


figure(3); 
subplot(1,3,3);
pwelch(outSig_PAPR_reduction, [], [], [], fs_baseband, 'centered');
title('PSD of Signal After PAPR Reduction'); 
ylim([-150,-80]); 
grid on;

outSig_c  = outSig(:);
outSig_clipped_c = outSig_clipped(:);
outSig_PAPR_reduction_c = outSig_PAPR_reduction(:);

PAPR_before    = 10*log10( max(abs(outSig_c)).^2  / mean(abs(outSig_c)).^2  );
PAPR_after     = 10*log10( max(abs(outSig_clipped_c)).^2 / mean(abs(outSig_clipped_c)).^2 );
PAPR_reduction = 10*log10( max(abs(outSig_PAPR_reduction_c)).^2 / mean(abs(outSig_PAPR_reduction_c)).^2 );


fprintf('PAPR Before Clipping: %.2f dB\n', PAPR_before);
fprintf('PAPR After Clipping: %.2f dB\n', PAPR_after);
fprintf('PAPR After PAPR Reduction: %.2f dB\n', PAPR_reduction);

%% [CCDF For All Stages]

signal_list = {outSig, outSig_clipped, outSig_PAPR_reduction};
labels = {'Before Clipping', 'After Clipping', 'After PAPR Reduction'};
legendEntries = {};

% [Plot] [6] CCDF 
figure(6);
for i = 1:length(signal_list)

    sig = signal_list{i};

    pm = powermeter(ComputeCCDF=true, PowerRange=50); 
    
    averagePower = pm(sig);  
    
    [relpower,prob] = ccdf(pm);
    
    plotCCDF(pm,GaussianReference=true)
    hold on
    xlabel('Relative Power (dB)');
    ylabel('Probability (%)');
  
    grid on;
    title('CCDF Measurement');
    legendEntries{end+1} = labels{i};
end

ylim([1e-3 1e2]);
xlim([-5 12]);
xticks(-5:0.5:12);

allLines = findobj(gca, 'Type', 'line');
ccdfLines = allLines(strcmp({allLines.LineStyle}, '-'));    
gaussLines = allLines(strcmp({allLines.LineStyle}, '--'));  
gaussLine = gaussLines(1); 
legend([flipud(ccdfLines); gaussLine], [legendEntries, {'Gaussian reference'}]);

%% [Frequency up-conversion]
baseband_signal = outSig_PAPR_reduction;
% upsample
UP_baseband_signal = resample(baseband_signal, FS_baseband, fs_baseband);
fc = 1e9; % Carrier frequency of 1 GHz
t = (0:length(UP_baseband_signal)-1).' / FS_baseband; % Time vector for the waveform
upconverted_outSig = real(UP_baseband_signal .* exp(1i * 2 * pi * fc * t));
figure(7);
pwelch(upconverted_outSig, [], [], [], FS_baseband, 'centered');

%% MZM
% E/O Conversion: MZM
Vpi = 4;        % MZM Vpi
Pavg = 1;        % Average optical power 

Vin = upconverted_outSig;
Vin_scaled = Vin / max(abs(Vin)) * (0.4*Vpi);

Popt = Pavg/2 * (1 + cos(pi * (Vin_scaled + Vpi/2) / Vpi));
% Optical Fibre Loss 
loss_dB     = 0.5;                               
loss_factor = 10^(-loss_dB / 10);
Popt = Popt * loss_factor;

% O/E Conversion
responsivity = 0.8;                              
electrical_signal = responsivity * Popt;         

% mm-wave Power Amplifier
gain_dB = 40;                                 
gain = 10^(gain_dB / 10);
amplified_signal = electrical_signal * sqrt(gain);

%% RF Downconversion
t2 = (0:numel(amplified_signal)-1).' / FS_baseband;
downconverted_outSig = amplified_signal.*exp(-1i * 2 * pi * fc * t2);
figure(8);
pwelch(downconverted_outSig, [], [], [], FS_baseband, 'centered');

%% Baseband Low-Pass Filtering
HdN02 = filterN02;
downconverted_outSig_filtered = filter(HdN02, downconverted_outSig);

L = round(FS_baseband / fs_baseband);  

gdN02 = round(mean(grpdelay(HdN02,1024)));      
TOTAL_DELAY = gdN01*L + gdN02;
downconverted_outSig_filtered = [downconverted_outSig_filtered(TOTAL_DELAY+1:end); ...
                                 zeros(TOTAL_DELAY,1)];
figure(9);
pwelch(downconverted_outSig_filtered, [], [], [], FS_baseband, 'centered');

%% DC remove
downconverted_outSig_filtered = downconverted_outSig_filtered - mean(downconverted_outSig_filtered);

%% resample
re_signal = resample(downconverted_outSig_filtered, fs_baseband, FS_baseband);
%re_signal = UP_baseband_signal(3:50:end);
%re_signal = downconverted_outSig_filtered(1:50:end);
figure(10)
pwelch(re_signal,[],[],[],'centered')
title('Basedband spectrum after down-sampling')
figure(11)
plot([real(re_signal) imag(re_signal)])
title('I and Q of the down-sampled time-domain signal')
%% OFDM Demodulation

data_out = ofdmdemod(re_signal, nfft, cplen, 0, nullIdx);
% Constellation Diagram
const_xlim = 200;

constel = comm.ConstellationDiagram( ...
    'ColorFading', true, ...
    'ShowTrajectory', false, ...
    'ShowReferenceConstellation', false, ...
    'XLimits', [-const_xlim const_xlim], ...
    'YLimits', [-const_xlim const_xlim]);

constel(data_out(:));
release(constel);

%% [Demodulation]

% EVM of outSig(original OFDM signal)
qamSym_origin = qamSym(:);

outSig_1 = ofdmdemod(outSig, nfft, cplen, 0, nullIdx);
outSig_2 = outSig_1(:);
evm_outSig = sqrt(mean(abs(outSig_2 - qamSym_origin).^2)) / sqrt(mean(abs(qamSym_origin).^2)) * 100;
fprintf('EVM of outSig: %.6f %%\n', evm_outSig);

% EVM of outSig_clipped(OFDM signal after clipping)
outSig_clipped_1 = ofdmdemod(outSig_clipped, nfft, cplen, 0, nullIdx);
outSig_clipped_2 = outSig_clipped_1(:);
evm_outSig_clipped = sqrt(mean(abs(outSig_clipped_2 - qamSym_origin).^2)) / sqrt(mean(abs(qamSym_origin).^2)) * 100;
fprintf('EVM of outSig_clipped: %.6f %%\n', evm_outSig_clipped);

% EVM of outSig_PAPR_reduction(final OFDM signal for transmission)
delay = 59;  
outSig_PAPR_reduction_sync = outSig_PAPR_reduction(delay+1:end);    

N = floor(length(outSig_PAPR_reduction_sync) / nfft) * nfft;
outSig_PAPR_reduction_sync = outSig_PAPR_reduction_sync(1:N);

outSig_PAPR_reduction_1 = ofdmdemod(outSig_PAPR_reduction_sync, nfft, cplen, 0, nullIdx);
outSig_PAPR_reduction_2 = outSig_PAPR_reduction_1(:);

qamSym_ref = qamSym_origin(1:length(outSig_PAPR_reduction_2));

evm_outSig_PAPR_reduction = sqrt(mean(abs(outSig_PAPR_reduction_2 - qamSym_ref).^2)) / ...
                            sqrt(mean(abs(qamSym_ref).^2)) * 100;

fprintf('EVM of outSig_PAPR_reduction: %.6f %%\n', evm_outSig_PAPR_reduction);

rx = data_out(:);                
tx = qamSym_origin(1:length(rx));

alpha = (tx' * rx) / (tx' * tx); 
rx_eq = rx / alpha;              

evm_final_eq = sqrt(mean(abs(rx_eq - tx).^2)) / sqrt(mean(abs(tx).^2)) * 100;
fprintf('EVM of final (after scalar EQ): %.3f %%\n', evm_final_eq);

%% [Constellation Diagram]

scatterplot(outSig_2(1:1000));
title('outSig');

scatterplot(outSig_clipped_2(1:1000));
title('outSig clipped');

scatterplot(outSig_PAPR_reduction_2(1:1000));
title('outSig PAPR reduction');

scatterplot(data_out(1:1000));
title('final');

data_out_c = data_out(:);

PAPR_final = 10*log10( max(abs(re_signal)).^2 / mean(abs(re_signal)).^2 );
fprintf('PAPR final: %.2f dB\n', PAPR_final);