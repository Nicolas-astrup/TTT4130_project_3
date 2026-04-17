%% Task 4 – Bit-Interleaved Coded Modulation (BICM)                                                                                                          
  clear; clc;                                                                                                                                                  
                                                                                                                                                               
  %--- Parameters ---
  msg_len   = 100;
  trellis   = poly2trellis(3, [7 5]);                                                                                                                          
  M         = 16;         % 16-QAM → 4 bits per symbol
  bps       = log2(M);    % = 4                                                                                                                                
  code_rate = 0.5;        % rate-1/2 convolutional code                                                                                                        
  Nrows     = 8;          % interleaver rows → minimum separation = Nrows
  Ncols     = 26;         % interleaver columns  (8×26 = 208 coded bits)                                                                                       
  tblen     = 35;                                                                                                                                            
                                                                                                                                                               
  crcGen = comm.CRCGenerator([1 0 0 1 1]);   % adds 4 CRC bits → 104 bits total                                                                                
  crcDet = comm.CRCDetector( [1 0 0 1 1]);
                                                                                                                                                               
  fprintf('Interleaver: %d rows × %d cols  →  min. separation = %d bits\n\n', ...                                                                              
      Nrows, Ncols, Nrows);
                                                                                                                                                               
  %--- Evaluate at 5 dB and 10 dB Eb/N0 ---                                                                                                                    
  for snr_dB = [5, 10]
                                                                                                                                                               
      %-- 1. Message + CRC                                                                                                                                     
      msg     = randi([0 1], msg_len, 1);
      msg_crc = crcGen(msg);                  % 104×1                                                                                                          
                                                                                                                                                               
      %-- 2. Convolutional encoding  (rate 1/2 → 208 bits)
      coded = convenc(msg_crc, trellis);      % 208×1                                                                                                          
                                                                                                                                                             
      %-- 3. Interleaving                                                                                                                                      
      %   Writes bits into 8×26 matrix row-by-row, reads out column-by-column.
      %   Adjacent encoded bits are separated by Nrows = 8 positions after interleaving.                                                                       
      coded_intrl = matintrlv(coded, Nrows, Ncols);                                                                                                            
                                                                                                                                                               
      %-- 4. 16-QAM modulation  (every 4 bits → 1 complex symbol → 52 symbols)                                                                                 
      tx = qammod(coded_intrl, M, 'InputType', 'bit', 'UnitAveragePower', true);                                                                             
                                                                                                                                                               
      %-- 5. AWGN channel                                                                                                                                    
      %   Es/N0 = Eb/N0 + 10·log10(bits_per_symbol × code_rate)                                                                                                
      %         = Eb/N0 + 10·log10(4 × 0.5) = Eb/N0 + 3 dB                                                                                                   
      EsN0_dB = snr_dB + 10*log10(bps * code_rate);                                                                                                            
      rx = awgn(tx, EsN0_dB);                                                                                                                                  
                                                                                                                                                               
      %-- 6. 16-QAM demodulation (hard decision)                                                                                                               
      rx_bits = qamdemod(rx, M, 'OutputType', 'bit', 'UnitAveragePower', true);                                                                              
                                                                                                                                                               
      %-- 7. De-interleaving                                                                                                                                 
      rx_deintrl = matdeintrlv(rx_bits, Nrows, Ncols);
                                                                                                                                                               
      %-- 8. Viterbi decoding (hard decision)
      decoded = vitdec(rx_deintrl, trellis, tblen, 'trunc', 'hard');                                                                                           
      decoded = decoded(1:length(msg_crc));   % 104 bits                                                                                                       
  
      %-- 9. CRC check                                                                                                                                         
      [recovered_msg, crc_error] = crcDet(decoded);                                                                                                          
                                                                                                                                                               
      %-- 10. BER
      bit_errors = sum(recovered_msg ~= msg);                                                                                                                  
      ber        = bit_errors / msg_len;                                                                                                                     
                                                                                                                                                               
      if crc_error, crc_str = 'FAIL'; else, crc_str = 'PASS'; end
                                                                                                                                                               
      fprintf('Eb/N0 = %2d dB | BER = %.4f  (%d/%d errors) | CRC: %s\n', ...                                                                                   
          snr_dB, ber, bit_errors, msg_len, crc_str);
  end                          