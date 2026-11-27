clc; clear; close all;
fprintf('ЛАБОРАТОРНАЯ РАБОТА №4 \n\n');

%% ПУНКТ 1. Исследование зависимости CCDF пик-фактора от кол-ва поднесущих

fprintf('[Пункт 1/5] Исследование зависимости CCDF от N_sub...\n');

N_sub_list_p1 = [64, 128, 256, 512, 1024, 2048]; 
Nfft_p1 = 2048;
M_inf_p1 = 16;
N_OFDM_sym_p1 = 20000; 

x_grid_p1 = 0:0.1:12; 
f_ccdf_grid_p1 = nan(length(x_grid_p1), length(N_sub_list_p1));

for m = 1:length(N_sub_list_p1)
    curr_N_sub = N_sub_list_p1(m);

    data_idx = randi([0 M_inf_p1-1], curr_N_sub, N_OFDM_sym_p1);
    mod_data = qammod(data_idx, M_inf_p1, 'UnitAveragePower', true);
    
    frame_f = zeros(Nfft_p1, N_OFDM_sym_p1);
    start_idx = Nfft_p1/2 - curr_N_sub/2 + 1;
    idx_range = start_idx : (start_idx + curr_N_sub - 1);
    frame_f(idx_range, :) = mod_data;
    
    tx_signal = ifft(ifftshift(frame_f, 1), Nfft_p1, 1);
    sig_vector = tx_signal(:);
    
    pwr = abs(sig_vector).^2;
    p_mean = mean(pwr);
    inst_papr_dB = 10 * log10(pwr / p_mean);
    
    [f_cdf, x_val] = ecdf(inst_papr_dB);
    f_ccdf_raw = 1 - f_cdf;
    
    valid_mask = isfinite(x_val);
    x_val = x_val(valid_mask);
    f_ccdf_raw = f_ccdf_raw(valid_mask);
    
    [xu, ~, ic] = unique(x_val);
    fu_max = accumarray(ic, f_ccdf_raw, [], @max);
    
    ccdf_curve = interp1(xu, fu_max, x_grid_p1, 'linear', NaN);
    f_ccdf_grid_p1(:, m) = ccdf_curve;
end

figure('Color', [0.94 0.94 0.94], 'Name', 'Пункт 1: CCDF от N_sub');
hold on; box on; 
colors_p1 = lines(length(N_sub_list_p1));
marker_step = 8; 
plot_handles_p1 = [];

for m = 1:length(N_sub_list_p1)
    y_data = f_ccdf_grid_p1(:, m);
    x_data = x_grid_p1;
    mask = ~isnan(y_data) & (y_data > 0);
    
    h = semilogy(x_data(mask), y_data(mask), 'LineWidth', 1.5, 'Color', colors_p1(m,:));
    plot_handles_p1 = [plot_handles_p1; h]; 
    
    valid_indices = find(mask);
    marker_indices = valid_indices(1:marker_step:end);
    semilogy(x_data(marker_indices), y_data(marker_indices), ...
        'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 5, ...
        'MarkerEdgeColor', colors_p1(m,:), 'MarkerFaceColor', 'none');
end
hold off;
set(gca, 'YScale', 'log', 'FontSize', 10, 'Box', 'on', 'Color', 'w');
xlabel('PAPR, дБ', 'FontSize', 11);
ylabel('CCDF', 'FontSize', 11);
title('Пункт 1: Зависимость CCDF пик-фактора от кол-ва ненулевых поднесущих', 'FontSize', 11);
xlim([0 12]); ylim([1e-8 1]);
legend_str_p1 = arrayfun(@(x) sprintf('N_{data} = %d', x), N_sub_list_p1, 'UniformOutput', false);
legend(plot_handles_p1, legend_str_p1, 'Location', 'northeast', 'FontSize', 9);


%% ПУНКТ 2. Исследование влияния ограничения (Clipping) на CCDF, PSD и BER

fprintf('[Пункт 2/5] Исследование влияния Clipping (CR = Inf, 8, 6, 4 дБ)...\n');

Nfft_p2 = 2048;
N_sub_p2 = 1024;           
M_p2 = 16;                
CP_p2 = 128;             
CR_dB_list_p2 = [Inf, 8, 6, 4]; 
EbNo_dB_p2 = 0:2:18;
N_sym_ber_p2 = 3000;    
N_sym_stat_p2 = 2000;      

bits_stat_p2 = randi([0 1], log2(M_p2)*N_sub_p2*N_sym_stat_p2, 1);
data_stat_p2 = qammod(bits_stat_p2, M_p2, 'InputType', 'bit', 'UnitAveragePower', true);
frame_stat_p2 = zeros(Nfft_p2, N_sym_stat_p2);
idx_active_p2 = Nfft_p2/2 - N_sub_p2/2 + 1 : Nfft_p2/2 + N_sub_p2/2;
frame_stat_p2(idx_active_p2, :) = reshape(data_stat_p2, N_sub_p2, N_sym_stat_p2);

tx_raw_p2 = ifft(ifftshift(frame_stat_p2, 1), Nfft_p2, 1);
tx_vec_p2 = tx_raw_p2(:);
sigma_p2 = std(tx_vec_p2);

bits_ber_p2 = randi([0 1], log2(M_p2)*N_sub_p2*N_sym_ber_p2, 1);
data_ber_p2 = qammod(bits_ber_p2, M_p2, 'InputType', 'bit', 'UnitAveragePower', true);
frame_ber_p2 = zeros(Nfft_p2, N_sym_ber_p2);
frame_ber_p2(idx_active_p2, :) = reshape(data_ber_p2, N_sub_p2, N_sym_ber_p2);
tx_ber_raw_p2 = ifft(ifftshift(frame_ber_p2, 1), Nfft_p2, 1);

results_ccdf_p2 = cell(length(CR_dB_list_p2), 1);
results_psd_p2 = cell(length(CR_dB_list_p2), 1);
results_ber_p2 = zeros(length(EbNo_dB_p2), length(CR_dB_list_p2));
x_grid_ccdf_p2 = 0:0.1:14; 

for i = 1:length(CR_dB_list_p2)
    cr_val = CR_dB_list_p2(i);
    if isinf(cr_val)
        threshold = inf;
    else
        threshold = 10^(cr_val/20) * sigma_p2;
    end

    sig_proc = tx_vec_p2;
    if ~isinf(threshold)
        mask = abs(sig_proc) > threshold;       
        sig_proc(mask) = threshold .* (sig_proc(mask) ./ abs(sig_proc(mask)));
    end
    
    pwr_inst = abs(sig_proc).^2;
    p_mean = mean(pwr_inst);
    papr_inst_dB = 10 * log10(pwr_inst / p_mean);
    
    [f_cdf, x_val] = ecdf(papr_inst_dB);
    f_ccdf = 1 - f_cdf;
    
    % Очистка от NaN и Inf
    valid_mask = isfinite(x_val) & isfinite(f_ccdf);
    x_val = x_val(valid_mask);
    f_ccdf = f_ccdf(valid_mask);
    
    [xu, unique_idx] = unique(x_val);
    fu = f_ccdf(unique_idx);
    results_ccdf_p2{i} = interp1(xu, fu, x_grid_ccdf_p2, 'pchip', NaN);
    
    [pxx, f_freq] = pwelch(sig_proc, hamming(2048), 1024, 2048, 1, 'centered');
    results_psd_p2{i} = [f_freq, 10*log10(pxx)];
    
    sig_ber_loop = tx_ber_raw_p2;
    if ~isinf(threshold)
        mask_ber = abs(sig_ber_loop) > threshold;
        sig_ber_loop(mask_ber) = threshold .* (sig_ber_loop(mask_ber) ./ abs(sig_ber_loop(mask_ber)));
    end
    
    sig_ber_cp = [sig_ber_loop(end-CP_p2+1:end, :); sig_ber_loop];
    
    for k = 1:length(EbNo_dB_p2)
        snr_val = EbNo_dB_p2(k) + 10*log10(log2(M_p2)) + 10*log10(N_sub_p2/Nfft_p2);
        rx_sig = awgn(sig_ber_cp, snr_val, 'measured');
        rx_no_cp = rx_sig(CP_p2+1:end, :);
        rx_fft = fftshift(fft(rx_no_cp, Nfft_p2, 1), 1);
        rx_syms = rx_fft(idx_active_p2, :);
        rx_bits = qamdemod(rx_syms(:), M_p2, 'OutputType', 'bit', 'UnitAveragePower', true);
        [~, results_ber_p2(k, i)] = biterr(bits_ber_p2, rx_bits);
    end
end

colors_p2 = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.494 0.184 0.556];
markers_p2 = {'none', 'o', 's', 'd'};
line_styles_p2 = {'-', '-', '-', '-'};

figure('Color', [0.94 0.94 0.94], 'Name', 'Пункт 2: Спектр (PSD)');
hold on; box on; 
for i = 1:length(CR_dB_list_p2)
    plot(results_psd_p2{i}(:,1), results_psd_p2{i}(:,2), ...
        'LineWidth', 1.5, 'Color', colors_p2(i,:), 'DisplayName', get_leg(CR_dB_list_p2(i)));
end
xlabel('Нормированная частота'); ylabel('PSD (дБ/Гц)');
title('Пункт 2: Энергетический спектр при ограничении (PSD)');
ylim([-100 -10]); xlim([-0.5 0.5]);
legend('show', 'Location', 'south'); set(gca, 'Color', 'w');

figure('Color', [0.94 0.94 0.94], 'Name', 'Пункт 2: CCDF мощности');
hold on; box on; set(gca, 'YScale', 'log'); 
plot_handles_p2 = [];
for i = 1:length(CR_dB_list_p2)
    y_data = results_ccdf_p2{i};
    mask_plot = ~isnan(y_data) & (y_data > 0);
    h = semilogy(x_grid_ccdf_p2(mask_plot), y_data(mask_plot), ...
        'LineWidth', 2, 'Color', colors_p2(i,:), 'LineStyle', line_styles_p2{i});
    plot_handles_p2 = [plot_handles_p2; h]; 
    if ~strcmp(markers_p2{i}, 'none')
        idx_mark = 1:10:length(x_grid_ccdf_p2); 
        idx_mark = idx_mark(mask_plot(idx_mark)); 
        semilogy(x_grid_ccdf_p2(idx_mark), y_data(idx_mark), 'LineStyle', 'none', ...
            'Marker', markers_p2{i}, 'MarkerEdgeColor', colors_p2(i,:), 'MarkerSize', 6);
    end
end
xlabel('PAPR, дБ'); ylabel('CCDF');
title('Пункт 2: CCDF при разном уровне ограничения');
xlim([0 12]); ylim([1e-5 1]);
legend(plot_handles_p2, arrayfun(@(x) get_leg(x), CR_dB_list_p2, 'UniformOutput', false), 'Location', 'southwest');
set(gca, 'Color', 'w');

figure('Color', [0.94 0.94 0.94], 'Name', 'Пункт 2: Помехоустойчивость (BER)');
semilogy(EbNo_dB_p2, berawgn(EbNo_dB_p2, 'qam', M_p2), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Theory');
hold on; box on; 
for i = 1:length(CR_dB_list_p2)
    semilogy(EbNo_dB_p2, results_ber_p2(:, i), 'LineWidth', 1.5, ...
        'Color', colors_p2(i,:), 'Marker', markers_p2{i}, 'MarkerSize', 6, 'DisplayName', get_leg(CR_dB_list_p2(i)));
end
xlabel('E_b / N_0 (дБ)'); ylabel('BER');
title('Пункт 2: Помехоустойчивость (BER) с учетом ограничений');
ylim([1e-5 1]); xlim([0 18]);
legend('show', 'Location', 'southwest'); set(gca, 'Color', 'w');


%% ПУНКТ 3. Итерационное ограничение и фильтрация (Clipping and Filtering - C&F)

fprintf('[Пункт 3/5] Итерационное ограничение и фильтрация (C&F, CR = 4 дБ)...\n');

Nfft_p3 = 2048;
N_sub_p3 = 1024;           
M_p3 = 16;                 
CP_p3 = 128;               
CR_dB_p3 = 4;              
Iterations_list_p3 = [1, 3, 5]; 
Scenarios_p3 = [0, Iterations_list_p3];

x_grid_ccdf_p3 = 0:0.1:12; 
EbNo_dB_p3 = 0:2:12;       
N_sym_stat_p3 = 2000;      
N_sym_ber_p3 = 3000;

start_idx_p3 = Nfft_p3/2 - N_sub_p3/2 + 1;
idx_active_p3 = start_idx_p3 : (start_idx_p3 + N_sub_p3 - 1);
filter_mask_p3 = false(Nfft_p3, 1);
filter_mask_p3(idx_active_p3) = true;

bits_stat_p3 = randi([0 1], log2(M_p3)*N_sub_p3*N_sym_stat_p3, 1);
data_stat_p3 = qammod(bits_stat_p3, M_p3, 'InputType', 'bit', 'UnitAveragePower', true);
frame_stat_p3 = zeros(Nfft_p3, N_sym_stat_p3);
frame_stat_p3(idx_active_p3, :) = reshape(data_stat_p3, N_sub_p3, N_sym_stat_p3);
tx_stat_raw_p3 = ifft(ifftshift(frame_stat_p3, 1), Nfft_p3, 1);

sigma_p3 = std(tx_stat_raw_p3(:));
threshold_p3 = 10^(CR_dB_p3/20) * sigma_p3;

bits_ber_p3 = randi([0 1], log2(M_p3)*N_sub_p3*N_sym_ber_p3, 1);
data_ber_p3 = qammod(bits_ber_p3, M_p3, 'InputType', 'bit', 'UnitAveragePower', true);
frame_ber_p3 = zeros(Nfft_p3, N_sym_ber_p3);
frame_ber_p3(idx_active_p3, :) = reshape(data_ber_p3, N_sub_p3, N_sym_ber_p3);
tx_ber_raw_p3 = ifft(ifftshift(frame_ber_p3, 1), Nfft_p3, 1);

results_ccdf_p3 = {};
results_psd_p3 = {};
results_ber_p3 = zeros(length(EbNo_dB_p3), length(Scenarios_p3));

for i = 1:length(Scenarios_p3)
    n_iter = Scenarios_p3(i);
    if n_iter == 0
        tx_proc = tx_stat_raw_p3;
        tx_ber_proc = tx_ber_raw_p3;
    else
        current_sig = tx_stat_raw_p3;
        for k = 1:n_iter
            mask = abs(current_sig) > threshold_p3;
            current_sig(mask) = threshold_p3 .* (current_sig(mask) ./ abs(current_sig(mask)));
            spec = fftshift(fft(current_sig, Nfft_p3, 1), 1);
            spec(~filter_mask_p3, :) = 0; 
            current_sig = ifft(ifftshift(spec, 1), Nfft_p3, 1);
        end
        tx_proc = current_sig;
        
        current_ber = tx_ber_raw_p3;
        for k = 1:n_iter
            mask = abs(current_ber) > threshold_p3;
            current_ber(mask) = threshold_p3 .* (current_ber(mask) ./ abs(current_ber(mask)));
            spec = fftshift(fft(current_ber, Nfft_p3, 1), 1);
            spec(~filter_mask_p3, :) = 0;
            current_ber = ifft(ifftshift(spec, 1), Nfft_p3, 1);
        end
        tx_ber_proc = current_ber;
    end
    
    sig_vec = tx_proc(:);
    pwr = abs(sig_vec).^2;
    papr_inst = 10*log10(pwr/mean(pwr));
    
    [f_cdf, x_val] = ecdf(papr_inst);
    [xu, idx_u] = unique(x_val);
    fu = 1 - f_cdf(idx_u);
    results_ccdf_p3{i} = interp1(xu, fu, x_grid_ccdf_p3, 'pchip', NaN);
    
    [pxx, f_freq] = pwelch(sig_vec, hamming(2048), 1024, 2048, 1, 'centered');
    results_psd_p3{i} = [f_freq, 10*log10(pxx)];
    
    sig_cp = [tx_ber_proc(end-CP_p3+1:end, :); tx_ber_proc];
    for k = 1:length(EbNo_dB_p3)
        snr_val = EbNo_dB_p3(k) + 10*log10(log2(M_p3)) + 10*log10(N_sub_p3/Nfft_p3);
        rx_sig = awgn(sig_cp, snr_val, 'measured');
        rx_no_cp = rx_sig(CP_p3+1:end, :);
        rx_fft = fftshift(fft(rx_no_cp, Nfft_p3, 1), 1);
        rx_syms = rx_fft(idx_active_p3, :);
        rx_bits = qamdemod(rx_syms(:), M_p3, 'OutputType', 'bit', 'UnitAveragePower', true);
        [~, results_ber_p3(k, i)] = biterr(bits_ber_p3, rx_bits);
    end
end

colors_p3 = [0 0 0; 0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188]; 
markers_p3 = {'none', 'o', 's', '^'};
line_styles_p3 = {'--', '-', '-', '-'};

figure('Color', [0.94 0.94 0.94], 'Name', 'Пункт 3: C&F Спектр (PSD)');
hold on; box on; 
for i = 1:length(Scenarios_p3)
    plot(results_psd_p3{i}(:,1), results_psd_p3{i}(:,2), ...
        'LineWidth', 1.5, 'Color', colors_p3(i,:), 'LineStyle', line_styles_p3{i}, ...
        'DisplayName', get_lbl(Scenarios_p3(i)));
end
xlabel('Норм. частота'); ylabel('PSD (дБ/Гц)');
title(['Пункт 3: Спектр при C&F (CR = ' num2str(CR_dB_p3) ' дБ)']);
legend('show', 'Location', 'south'); ylim([-100 -10]); xlim([-0.5 0.5]);
set(gca, 'Color', 'w');

figure('Color', [0.94 0.94 0.94], 'Name', 'Пункт 3: C&F CCDF');
hold on; box on; set(gca, 'YScale', 'log');
for i = 1:length(Scenarios_p3)
    y_data = results_ccdf_p3{i};
    mask = ~isnan(y_data) & (y_data > 0);
    semilogy(x_grid_ccdf_p3(mask), y_data(mask), ...
        'LineWidth', 2, 'Color', colors_p3(i,:), 'LineStyle', line_styles_p3{i}, ...
        'DisplayName', get_lbl(Scenarios_p3(i)));
    if i > 1
        idx = 1:10:length(x_grid_ccdf_p3);
        idx = intersect(idx, find(mask));
        semilogy(x_grid_ccdf_p3(idx), y_data(idx), 'LineStyle', 'none', ...
            'Marker', markers_p3{i}, 'MarkerEdgeColor', colors_p3(i,:), 'MarkerSize', 6);
    end
end
xlabel('PAPR, дБ'); ylabel('CCDF');
title('Пункт 3: Снижение PAPR методом C&F по итерациям');
xlim([0 12]); ylim([1e-4 1]); legend('show', 'Location', 'southwest');
set(gca, 'Color', 'w');

figure('Color', [0.94 0.94 0.94], 'Name', 'Пункт 3: C&F BER');
semilogy(EbNo_dB_p3, berawgn(EbNo_dB_p3, 'qam', M_p3), 'k:', 'LineWidth', 1.5, 'DisplayName', 'Theory');
hold on; box on;
for i = 1:length(Scenarios_p3)
    semilogy(EbNo_dB_p3, results_ber_p3(:, i), ...
        'LineWidth', 1.5, 'Color', colors_p3(i,:), 'LineStyle', line_styles_p3{i}, ...
        'Marker', markers_p3{i}, 'DisplayName', get_lbl(Scenarios_p3(i)));
end
xlabel('Eb/N0, дБ'); ylabel('BER');
title(['Пункт 3: BER при C&F (CR = ' num2str(CR_dB_p3) ' дБ)']);
ylim([1e-5 1]); legend('show', 'Location', 'southwest'); set(gca, 'Color', 'w');


%% ПУНКТ 4. Метод SLM (Selected Mapping) с 4-фазным алфавитом

fprintf('[Пункт 4/5] Моделирование SLM метода (4 фазы)...\n');

Nfft_p4 = 2048;
N_sub_p4 = 1024;         
M_p4 = 16;                
U_list_p4 = [0, 4, 16, 128, 516]; 
N_sym_p4 = 1000;           
x_grid_p4 = 0:0.1:13;
results_ccdf_p4 = cell(length(U_list_p4), 1);
phase_alphabet_p4 = [1, 1j, -1, -1j]; 

bits_p4 = randi([0 1], log2(M_p4)*N_sub_p4*N_sym_p4, 1);
data_p4 = qammod(bits_p4, M_p4, 'InputType', 'bit', 'UnitAveragePower', true);
data_matrix_p4 = reshape(data_p4, N_sub_p4, N_sym_p4);
start_idx_p4 = Nfft_p4/2 - N_sub_p4/2 + 1;
idx_active_p4 = start_idx_p4 : (start_idx_p4 + N_sub_p4 - 1);

for k = 1:length(U_list_p4)
    U_req = U_list_p4(k);
    if U_req == 0, U_eff = 1; else, U_eff = U_req; end
    winner_signals = zeros(Nfft_p4 * N_sym_p4, 1);
    
    for s = 1:N_sym_p4
        X_curr = data_matrix_p4(:, s);
        if U_eff == 1
            P = ones(N_sub_p4, 1);
        else
            phase_indices = randi([1, 4], N_sub_p4, U_eff);
            P = phase_alphabet_p4(phase_indices);
            P(:, 1) = 1; 
        end
        X_candidates = repmat(X_curr, 1, U_eff) .* P;
        Frame_cands = zeros(Nfft_p4, U_eff);
        Frame_cands(idx_active_p4, :) = X_candidates;
        tx_cands = ifft(ifftshift(Frame_cands, 1), Nfft_p4, 1);
        
        pwr = abs(tx_cands).^2;
        papr_vals = max(pwr, [], 1) ./ mean(pwr, 1);
        [~, best_idx] = min(papr_vals);
        winner_signals((s-1)*Nfft_p4 + 1 : s*Nfft_p4) = tx_cands(:, best_idx);
    end
    
    pwr_inst = abs(winner_signals).^2;
    papr_dB = 10 * log10(pwr_inst / mean(pwr_inst));
    [f_cdf, x_val] = ecdf(papr_dB);
    f_ccdf_raw = 1 - f_cdf;
    valid = isfinite(x_val);
    [xu, idx_u] = unique(x_val(valid));
    results_ccdf_p4{k} = interp1(xu, f_ccdf_raw(valid(idx_u)), x_grid_p4, 'pchip', NaN);
end

figure('Color', [0.94 0.94 0.94], 'Name', 'Пункт 4: SLM Analysis');
hold on; box on; 
colors_p4 = [0 0 0; 1 0 0; 0 0 1; 0 0.5 0; 0.5 0 0.5]; 
markers_p4 = {'none', 'o', 's', '^', 'd'};
plot_handles_p4 = []; legend_str_p4 = {};

for k = 1:length(U_list_p4)
    y_data = results_ccdf_p4{k};
    mask = ~isnan(y_data) & (y_data > 0);
    if U_list_p4(k) == 0, lbl = 'Без SLM (Original)'; lw = 2.5; else, lbl = sprintf('SLM U=%d (4 фазы)', U_list_p4(k)); lw = 2.0; end
    legend_str_p4{end+1} = lbl; 
    h = semilogy(x_grid_p4(mask), y_data(mask), 'LineWidth', lw, 'Color', colors_p4(k,:));
    plot_handles_p4 = [plot_handles_p4; h]; 
    if k > 1
        idx = 1:8:length(x_grid_p4);
        idx = intersect(idx, find(mask));
        semilogy(x_grid_p4(idx), y_data(idx), 'LineStyle', 'none', ...
            'Marker', markers_p4{k}, 'MarkerEdgeColor', colors_p4(k,:), 'MarkerSize', 6);
    end
end
set(gca, 'YScale', 'log');
xlabel('PAPR, дБ'); ylabel('CCDF');
title('Пункт 4: Эффективность SLM метода с 4-фазным алфавитом');
xlim([0 12.5]); ylim([1e-5 1]);
legend(plot_handles_p4, legend_str_p4, 'Location', 'southwest'); set(gca, 'Color', 'w');


%% ПУНКТ 5. Метод PTS (Partial Transmit Sequence)

fprintf('[Пункт 5/5] Моделирование PTS метода...\n');

Nfft_p5 = 2048;
N_sub_p5 = 1024;          
M_p5 = 16;                 
V_list_p5 = [1, 2, 4, 8]; 
N_sym_p5 = 1000;           
x_grid_p5 = 0:0.1:13;
results_ccdf_p5 = cell(length(V_list_p5), 1);

bits_p5 = randi([0 1], log2(M_p5)*N_sub_p5*N_sym_p5, 1);
data_p5 = qammod(bits_p5, M_p5, 'InputType', 'bit', 'UnitAveragePower', true);
data_matrix_p5 = reshape(data_p5, N_sub_p5, N_sym_p5);
start_idx_p5 = Nfft_p5/2 - N_sub_p5/2 + 1;
all_active_idx_p5 = start_idx_p5 : (start_idx_p5 + N_sub_p5 - 1);

for k = 1:length(V_list_p5)
    V = V_list_p5(k);
    if V > 1
        num_combs = 2^(V-1);
        phase_combs = ones(num_combs, V);
        phase_combs(:, 2:end) = 1 - 2 * ff2n(V-1); 
    else
        phase_combs = 1;
    end
    
    winner_signals = zeros(Nfft_p5 * N_sym_p5, 1);
    
    for s = 1:N_sym_p5
        X_curr = data_matrix_p5(:, s);
        partial_time_signals = zeros(Nfft_p5, V);
        block_size = floor(N_sub_p5 / V);
        
        for v = 1:V
            if v < V, sub_idx_local = (v-1)*block_size + 1 : v*block_size;
            else, sub_idx_local = (v-1)*block_size + 1 : N_sub_p5; end
            
            frame_temp = zeros(Nfft_p5, 1);
            frame_temp(all_active_idx_p5(sub_idx_local)) = X_curr(sub_idx_local);
            partial_time_signals(:, v) = ifft(ifftshift(frame_temp, 1), Nfft_p5, 1);
        end
        
        min_papr = inf;
        best_sig = zeros(Nfft_p5, 1);
        for c = 1:size(phase_combs, 1)
            cand_sig = partial_time_signals * phase_combs(c, :).';
            pwr = abs(cand_sig).^2;
            curr_papr = max(pwr) / mean(pwr);
            if curr_papr < min_papr
                min_papr = curr_papr;
                best_sig = cand_sig;
            end
        end
        winner_signals((s-1)*Nfft_p5 + 1 : s*Nfft_p5) = best_sig;
    end
    
    pwr_inst = abs(winner_signals).^2;
    papr_dB = 10 * log10(pwr_inst / mean(pwr_inst));
    [f_cdf, x_val] = ecdf(papr_dB);
    f_ccdf = 1 - f_cdf;
    
    % Надежная фильтрация по бесконечным и NaN значениям
    valid = isfinite(x_val) & isfinite(f_ccdf);
    x_val = x_val(valid);
    f_ccdf = f_ccdf(valid);
    
    [xu, idx_u] = unique(x_val);
    fu = f_ccdf(idx_u);
    
    results_ccdf_p5{k} = interp1(xu, fu, x_grid_p5, 'pchip', NaN);
end

figure('Color', [0.94 0.94 0.94], 'Name', 'Пункт 5: PTS Analysis');
hold on; box on; 
colors_p5 = [0 0 0; 1 0 0; 0 0 1; 0 0.5 0]; 
markers_p5 = {'none', 'o', 's', '^'};
plot_handles_p5 = []; legend_str_p5 = {};

for k = 1:length(V_list_p5)
    y_data = results_ccdf_p5{k};
    mask = ~isnan(y_data) & (y_data > 0);
    if V_list_p5(k) == 1, lbl = 'Original (V=1)'; lw = 2.5; else, lbl = sprintf('PTS (V=%d)', V_list_p5(k)); lw = 2.0; end
    legend_str_p5{end+1} = lbl; 
    h = semilogy(x_grid_p5(mask), y_data(mask), 'LineWidth', lw, 'Color', colors_p5(k,:));
    plot_handles_p5 = [plot_handles_p5; h]; 
    if k > 1
        idx = 1:8:length(x_grid_p5);
        idx = intersect(idx, find(mask));
        semilogy(x_grid_p5(idx), y_data(idx), 'LineStyle', 'none', ...
            'Marker', markers_p5{k}, 'MarkerEdgeColor', colors_p5(k,:), 'MarkerSize', 6);
    end
end
set(gca, 'YScale', 'log');
xlabel('PAPR, дБ'); ylabel('CCDF');
title('Пункт 5: Эффективность PTS метода при разном числе блоков V');
xlim([0 12]); ylim([1e-5 1]);
legend(plot_handles_p5, legend_str_p5, 'Location', 'southwest'); set(gca, 'Color', 'w');

fprintf('\n ВСЕ РАСЧЕТЫ И ГРАФИКИ УСПЕШНО ВЫПОЛНЕНЫ! \n');


%% =========================================================================
% ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ 
% =========================================================================
function s = get_leg(val)
    if isinf(val)
        s = 'No Clipping';
    else
        s = sprintf('CR = %d dB', val);
    end
end

function s = get_lbl(n)
    if n == 0
        s = 'No Clip';
    else
        s = sprintf('%d Iterations', n);
    end
end