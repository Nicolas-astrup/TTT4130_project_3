%% Task 4 – Interleaved Coded 16-QAM over AWGN
clear; clc;

%% Parameters
msg_len    = 1000;          % message bits
crc_poly   = [1 0 0 1 1];  % CRC-4 polynomial
trellis    = poly2trellis(3, [7 5]);  % K=3, rate 1/2, octal [7,5]
tblen      = 35;            % Viterbi traceback length
M          = 16;            % 16-QAM
bps        = log2(M);       % bits per symbol = 4
interl_rows = 10;           % interleaver: rows (= burst-error separation)
SNR_dB_vec = 0:2:18;       % SNR sweep

%% CRC objects
crcGen = comm.CRCGenerator(crc_poly);
crcDet = comm.CRCDetector(crc_poly);

%% Pre-compute interleaver column count from encoded block size
n_crc      = msg_len + length(crc_poly) - 1;  % bits after CRC
n_enc      = 2 * n_crc;                        % bits after rate-1/2 conv enc
% Pad encoded length to a multiple of (interl_rows * bps) so QAM works cleanly
pad_to     = interl_rows * bps;
n_padded   = ceil(n_enc / pad_to) * pad_to;
interl_cols = n_padded / interl_rows;

fprintf('=== Task 4 System Parameters ===\n');
fprintf('Message bits    : %d\n', msg_len);
fprintf('After CRC       : %d bits\n', n_crc);
fprintf('After conv. enc : %d bits\n', n_enc);
fprintf('Padded to       : %d bits  (interleaver %dx%d)\n', ...
    n_padded, interl_rows, interl_cols);
fprintf('Min. burst separation (interleaver depth): %d bits\n\n', interl_rows);

%% SNR sweep
BER_uncoded  = zeros(size(SNR_dB_vec));
BER_coded    = zeros(size(SNR_dB_vec));

for idx = 1:length(SNR_dB_vec)
    EbN0_dB = SNR_dB_vec(idx);
    n_trials = 50;
    err_coded   = 0;
    total_bits  = 0;

    for trial = 1:n_trials
        %% --- TRANSMITTER ---
        % 1. Generate message
        msg = randi([0 1], msg_len, 1);

        % 2. CRC encode
        msg_crc = crcGen(msg);          % column vector, length n_crc

        % 3. Convolutional encode (flush registers)
        msg_flushed = [msg_crc; zeros(2, 1)];         % flush K-1=2 zeros
        enc = convenc(msg_flushed(:)', trellis);       % row vector, ~2*n_crc bits
        enc = enc(1:n_enc);                             % trim exactly

        % 4. Zero-pad to interleaver size
        enc_padded = [enc, zeros(1, n_padded - n_enc)];

        % 5. Block interleave: write row-by-row, read column-by-column
        %    → a burst of 'interl_rows' errors becomes spread by 'interl_cols'
        mat_tx     = reshape(enc_padded, interl_cols, interl_rows);
        interleaved = reshape(mat_tx', 1, []);          % read row-by-row after transpose

        % 6. 16-QAM modulation (Gray coded)
        bits_for_qam = reshape(interleaved, bps, [])'; % Nsym x 4 matrix
        syms = qammod(bi2de(bits_for_qam, 'left-msb'), M, ...
            'gray', 'InputType', 'integer', 'UnitAveragePower', true);

        %% --- CHANNEL ---
        % Eb/N0 → SNR per symbol: SNR_sym = Eb/N0 * bps * code_rate
        % code_rate = 1/2, so SNR_sym_dB = EbN0_dB + 10*log10(bps/2)
        SNR_sym_dB = EbN0_dB + 10*log10(bps * 0.5);
        rx_syms    = awgn(syms, SNR_sym_dB, 'measured');

        %% --- RECEIVER ---
        % 7. 16-QAM demodulate (hard decision)
        rx_int  = qamdemod(rx_syms, M, 'gray', ...
            'OutputType', 'integer', 'UnitAveragePower', true);
        rx_bits = de2bi(rx_int, bps, 'left-msb');
        rx_bits = reshape(rx_bits', 1, []);  % row vector

        % 8. De-interleave (inverse of step 5)
        mat_rx      = reshape(rx_bits, interl_rows, interl_cols);
        deinterl    = reshape(mat_rx', 1, []);          % read column-by-column after transpose
        deinterl_enc = deinterl(1:n_enc);               % strip padding

        % 9. Viterbi decode (hard decision)
        dec = vitdec(deinterl_enc, trellis, tblen, 'trunc', 'hard');
        dec = dec(1:n_crc)';                            % column vector

        % 10. CRC check
        [recovered, crc_err] = crcDet(dec);

        % BER: compare recovered message to original
        n_err = sum(recovered ~= msg);
        err_coded  = err_coded + n_err;
        total_bits = total_bits + msg_len;
    end

    BER_coded(idx) = err_coded / total_bits;

    % Uncoded reference: raw 16-QAM BER (theoretical)
    BER_uncoded(idx) = berawgn(EbN0_dB, 'qam', M);

    fprintf('SNR = %2d dB  |  BER coded = %.2e  |  BER uncoded = %.2e\n', ...
        EbN0_dB, BER_coded(idx), BER_uncoded(idx));
end

%% Plot
figure;
semilogy(SNR_dB_vec, BER_uncoded, 'b--o', 'LineWidth', 1.5, 'DisplayName', 'Uncoded 16-QAM (theory)');
hold on;
semilogy(SNR_dB_vec, max(BER_coded, 1e-6), 'r-s', 'LineWidth', 1.5, 'DisplayName', 'Coded: CRC + Conv + Interl + 16-QAM');
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('BER');
title('Task 4 – BER vs SNR: Coded 16-QAM');
legend('Location', 'southwest');
ylim([1e-6 1]);

%% Print results at requested SNR points
fprintf('\n--- Results at specific SNR points ---\n');
for snr_check = [5, 10]
    [~, i] = min(abs(SNR_dB_vec - snr_check));
    fprintf('At %2d dB: BER coded = %.4e,  BER uncoded = %.4e\n', ...
        snr_check, BER_coded(i), BER_uncoded(i));
end