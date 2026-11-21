%%
clear; close all; clc;

cfg = lab3_cfg();
if ~exist(cfg.figDir, 'dir'); mkdir(cfg.figDir); end

if strcmpi(getenv('LAB3_DEBUG'), '1')
    lab3_debug(cfg);
    return;
end

diary(fullfile(cfg.outDir, 'Lab3_results.txt'));
diary on;
tAll = tic;

fprintf('=============================================================\n');
fprintf(' ЛАБОРАТОРНАЯ РАБОТА №3: синхронизация OFDM в канале с АБГШ\n');
fprintf('=============================================================\n');
fprintf(' Nfft = %d, Ncp = %d (%.0f%% от Nfft), кадр: преамбула + %d символов\n', ...
    cfg.Nfft, cfg.Ncp, 100*cfg.Ncp/cfg.Nfft, cfg.nDataSym);
fprintf(' Данные: КАМ-16 на %d поднесущих, пилоты: КАМ-4 на %d поднесущих\n', ...
    numel(cfg.dataIdx), numel(cfg.pilotIdx));
fprintf(' Оценка BER по символам %d...%d, бит в кадре: %d\n', ...
    cfg.evalSym(1), cfg.evalSym(end), cfg.qamBits*numel(cfg.dataIdx)*numel(cfg.evalSym));
fprintf(' Eb/N0 = %s дБ, кадров на точку <= %d\n\n', mat2str(cfg.snrDb), cfg.nFrames);

S = struct('tag', {}, 'snr', {}, 'ber', {}, 'err', {}, 'bits', {});

%%  ЧАСТЬ 0. Демонстрация временного и частотного сдвига

fprintf('--- Часть 0. Временной и частотный сдвиг ---------------------\n');

rng(cfg.seed);
txFrame = lab_frame(cfg);

% (а) временной сдвиг: кадр запаздывает на delta отсчётов
deltaDemo = 300;
rxRef = lab_channel(cfg, txFrame, 0,          0, cfg.snrDemo);
rxDel = lab_channel(cfg, txFrame, deltaDemo,  0, cfg.snrDemo);

% (б,в) фазовая ошибка от временного сдвига внутри ЦП: до и после эквалайзера
%       созвездие показываем при высоком ОСШ, иначе точки "тонут" в шуме
deltaCP  = -cfg.Ncp/8;
snrConst = 22;                                  % Eb/N0 для демонстрации созвездий
rxCP = lab_channel(cfg, txFrame, deltaCP, 0, snrConst);
win0 = cfg.padLen + 1 + cfg.Ncp + (cfg.evalSym(1))*cfg.Nsym;
Seg  = rxCP(win0 + (0:cfg.Nfft-1).');
Ycp  = fft(Seg)/sqrt(cfg.Nfft);
Hcp  = lab_channel_estimate(cfg, Ycp);
YcpEq = Ycp ./ Hcp;

h = new_figure(cfg, 'Сдвиги сигнала', 'Временной и частотный сдвиг OFDM-сигнала');
tiledlayout(h, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% (а) временной сдвиг виден как смещение максимума корреляции по ЦП
[~, mRef, axMet] = cp_metric(rxRef, cfg, 1, cfg.padLen + 1);
[~, mDel]        = cp_metric(rxDel, cfg, 1, cfg.padLen + 1);
nexttile;
plot(axMet, mRef, 'Color', [0 0.45 0.74], 'LineWidth', 1.2, 'DisplayName', '\delta = 0');
hold on; grid on;
plot(axMet, mDel, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2, ...
     'DisplayName', sprintf('\\delta = +%d', deltaDemo));
[~, iRef] = max(mRef);  [~, iDel] = max(mDel);
plot(axMet(iRef), mRef(iRef), 'o', 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74], ...
     'MarkerSize', 7, 'HandleVisibility', 'off');
plot(axMet(iDel), mDel(iDel), 'o', 'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10], ...
     'MarkerSize', 7, 'HandleVisibility', 'off');
xlabel('d, отсчёты'); ylabel('M(d)');
title(sprintf('(а) \\delta = +%d: пик смещён на %d', deltaDemo, ...
      round(axMet(iDel) - axMet(iRef))));
legend('Location', 'northeast'); style_axes();

nexttile;
plot(Ycp(cfg.dataIdx), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 3);
grid on;
xlim([-2 2]); ylim([-2 2]); axis square;
xlabel('I'); ylabel('Q');
title('(б) без эквалайзера');
style_axes();

nexttile;
plot(YcpEq(cfg.dataIdx), '.', 'Color', [0.47 0.67 0.19], 'MarkerSize', 3);
grid on; xlim([-2 2]); ylim([-2 2]); axis square;
xlabel('I'); ylabel('Q');
title('(в) после эквалайзера');
style_axes();
save_figure(cfg, h, 'fig_00_shifts');

fprintf('  Рис. 0: сдвиг кадра на %d отсчётов; сдвиг на %d отсчётов внутри ЦП\n', ...
    deltaDemo, deltaCP);
fprintf('         даёт только фазовый наклон по поднесущим (снимается эквалайзером).\n\n');

%% ЧАСТЬ 1. Временная синхронизация
fprintf('--- Часть 1. Временная синхронизация -------------------------\n');

% 1.1 Помехоустойчивость при временной расстройке (с эквалайзером и без)
%     delta = 0 | -Ncp/2 (внутри ЦП) | -(Ncp+64) (за пределами ЦП) | +96
deltaList = [0, -cfg.Ncp/8, -(cfg.Ncp + 64), 96];
deltaName = {'\delta = 0', sprintf('\\delta = -%d (|\\delta| < N_{CP})', cfg.Ncp/8), ...
             sprintf('\\delta = -%d (|\\delta| > N_{CP})', cfg.Ncp + 64), '\delta = +96'};
tags11 = cell(1, 2*numel(deltaList));
for iD = 1:numel(deltaList)
    for iEq = 1:2                       % 1 - с эквалайзером, 2 - без
        scen = scen_set(make_scen(), 'delta', deltaList(iD), ...
                        'useEq', (iEq == 1), 'timeSync', 'none');
        tag = sprintf('f11_d%d_eq%d', iD, iEq == 1);
        out = ber_sweep(cfg, scen, cfg.snrDb, sprintf('%s, %s', deltaName{iD}, ...
              ternary(iEq == 1, 'с экв.', 'без экв.')));
        S = add_curve(S, tag, out);
        tags11{(iD-1)*2 + iEq} = tag;
    end
end

h = new_figure(cfg, 'BER при временной расстройке', 'Помехоустойчивость при временной расстройке');
ax = axes(h); hold(ax, 'on'); grid(ax, 'on');
cols = [0 0.45 0.74; 0.85 0.33 0.10; 0.49 0.18 0.56; 0.47 0.67 0.19];
for iD = 1:numel(deltaList)
    for iEq = 1:2
        c = get_curve(S, tags11{(iD-1)*2 + iEq});
        st = struct('col', cols(iD,:), 'mk', 'os^dv', ...
                    'ls', ternary(iEq == 1, '-', '--'), 'open', iEq == 2);
        st.mk = st.mk(iD);
        plot_ber_curve(ax, c, st, sprintf('%s, %s', deltaName{iD}, ...
            ternary(iEq == 1, 'с эквалайзером', 'без эквалайзера')));
    end
end
plot_theory(ax, cfg);
finish_ber_axes(ax, cfg, sprintf('N_{CP} = %d, \\epsilon = 0', cfg.Ncp));
save_figure(cfg, h, 'fig_1_1_delta');

% 1.1б Помехоустойчивость в зависимости от величины временного сдвига
deltaSweep = [-512 -448 -384 -320 -256 -224 -192 -160 -128 -96 -64 -32 -16 0 ...
              16 32 64 96 128 192 256 320 384 448 512];
berSweep = zeros(numel(deltaSweep), 2);
for iEq = 1:2
    for iD = 1:numel(deltaSweep)
        scen = scen_set(make_scen(), 'delta', deltaSweep(iD), ...
                        'useEq', (iEq == 1), 'timeSync', 'none');
        berSweep(iD, iEq) = ber_point(cfg, scen, cfg.snrDelta, cfg.nFramesDelta, 600);
    end
end
fprintf('  Развёртка по \\delta при Eb/N0 = %g дБ завершена\n', cfg.snrDelta);

h = new_figure(cfg, 'BER от временного сдвига', 'Помехоустойчивость от величины временного сдвига');
ax = axes(h); hold(ax, 'on'); grid(ax, 'on');
semilogy(ax, deltaSweep, berSweep(:,1), '-o', 'Color', [0 0.45 0.74], 'LineWidth', 1.5, ...
         'MarkerSize', 5, 'MarkerFaceColor', [0 0.45 0.74], 'DisplayName', 'с эквалайзером');
semilogy(ax, deltaSweep, berSweep(:,2), '--s', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, ...
         'MarkerSize', 5, 'MarkerFaceColor', 'w', 'DisplayName', 'без эквалайзера');
set(ax, 'YScale', 'log', 'XScale', 'linear');
ylim(ax, [1e-5 1]);
yl = ylim(ax);
plot(ax, [-cfg.Ncp -cfg.Ncp], yl, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot(ax, [0 0], yl, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
text(ax, -cfg.Ncp, yl(1)*8, ' -N_{CP}', 'FontSize', 10);
text(ax, 8, yl(1)*8, '\delta = 0', 'FontSize', 10);
xlabel(ax, 'временной сдвиг \delta, отсчёты');
ylabel(ax, 'BER');
title(ax, sprintf('Eb/N0 = %g дБ; слева от -N_{CP} и справа от 0 - межсимвольная интерференция', cfg.snrDelta));
legend(ax, 'Location', 'northwest');
style_axes(ax);
save_figure(cfg, h, 'fig_1_1b_delta_sweep');

% 1.2 Влияние числа символов когерентного накопления M на метрику
%     Пик метрики от M не зависит, а шумовой пьедестал падает как 1/M,
%     поэтому ось ординат - логарифмическая, иначе разница не видна.
Mlist = [1 2 4 8];
h = new_figure(cfg, 'Метрика ЦП', 'Метрика временной синхронизации по ЦП');
ax = axes(h); hold(ax, 'on'); grid(ax, 'on');
mStyles = {'-', '--', '-.', ':'};
hleg = [];
for iM = 1:numel(Mlist)
    [mAvg, axD] = metric_shape(cfg, make_scen(), Mlist(iM), cfg.snrDemo, cfg.nAvgMetric);
    mAvg(mAvg <= 0) = NaN;
    pl = plot(ax, axD, mAvg, 'LineStyle', mStyles{iM}, 'Color', cols(iM,:), 'LineWidth', 1.5, ...
              'DisplayName', sprintf('M = %d', Mlist(iM)));
    hleg = [hleg pl];
    far = abs(axD) > 500;                       % уровень вне главного лепестка
    fprintf(['  M = %d: пик метрики %.3f, средний шумовой пьедестал вне пика ' ...
             '%.5f (отношение %.1f)\n'], Mlist(iM), max(mAvg), ...
            mean(mAvg(far), 'omitnan'), max(mAvg)/mean(mAvg(far), 'omitnan'));
end
set(ax, 'YScale', 'log');
ylim(ax, [1e-5 0.3]);
xlabel('отсчёты относительно истинного начала кадра');
ylabel('M(d), отн. ед.');
title(sprintf('Метрика по ЦП при Eb/N0 = %d дБ (усреднение по %d кадрам)', ...
      cfg.snrDemo, cfg.nAvgMetric));
legend(hleg, 'Location', 'southeast'); style_axes();
save_figure(cfg, h, 'fig_1_2_metric_M');

% 1.3 Влияние M на помехоустойчивость (синхронизация по ЦП, эквалайзер)
%     (а) штатный ЦП (Ncp = 256): корреляционный выигрыш велик, метрика
%         надёжна, поэтому BER от M практически не зависит;
%     (б) короткий ЦП (Ncp = 32): выигрыш корреляции всего 15 дБ, при низком
%         ОСШ метрика ошибается, и накопление заметно улучшает приём.
tags13 = cell(1, numel(Mlist));
for iM = 1:numel(Mlist)
    scen = scen_set(make_scen(), 'deltaMode', 'uniform', 'timeSync', 'cp', ...
                    'M', Mlist(iM), 'useEq', true);
    tag = sprintf('f13_M%d', Mlist(iM));
    out = ber_sweep(cfg, scen, cfg.snrDb, sprintf('ЦП, M = %d', Mlist(iM)));
    S = add_curve(S, tag, out); tags13{iM} = tag;
end
outIdeal1 = ber_sweep(cfg, scen_set(make_scen(), 'deltaMode', 'uniform', ...
                      'timeSync', 'genie'), cfg.snrDb, 'точное положение кадра');
S = add_curve(S, 'f13_ideal', outIdeal1);

cfgShort = cfg;
cfgShort.Ncp = 32;  cfgShort.Nsym = cfgShort.Nfft + cfgShort.Ncp;
tags13b = cell(1, numel(Mlist));
for iM = 1:numel(Mlist)
    scen = scen_set(make_scen(), 'deltaMode', 'uniform', 'timeSync', 'cp', ...
                    'M', Mlist(iM), 'useEq', true);
    tag = sprintf('f13b_M%d', Mlist(iM));
    out = ber_sweep(cfgShort, scen, cfg.snrLow, sprintf('N_{CP}=32, M = %d', Mlist(iM)));
    S = add_curve(S, tag, out); tags13b{iM} = tag;
end
outIdeal2 = ber_sweep(cfgShort, scen_set(make_scen(), 'deltaMode', 'uniform', ...
                      'timeSync', 'genie'), cfg.snrLow, 'точное положение кадра');
S = add_curve(S, 'f13b_ideal', outIdeal2);

h = new_figure(cfg, 'BER против M', 'Влияние когерентного накопления на помехоустойчивость');
tiledlayout(h, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hold on; grid on;
for iM = 1:numel(Mlist)
    st = struct('col', cols(iM,:), 'mk', 'os^d', 'ls', '-', 'open', false);
    st.mk = st.mk(iM);
    plot_ber_curve(gca, get_curve(S, tags13{iM}), st, sprintf('M = %d', Mlist(iM)));
end
plot_ber_curve(gca, get_curve(S, 'f13_ideal'), ...
    struct('col', [0 0 0], 'mk', 'none', 'ls', ':', 'open', true), 'точное положение кадра');
plot_theory(gca, cfg);
finish_ber_axes(gca, cfg, sprintf('(а) N_{CP} = %d: кривые совпадают', cfg.Ncp));

nexttile;
hold on; grid on;
for iM = 1:numel(Mlist)
    st = struct('col', cols(iM,:), 'mk', 'os^d', 'ls', '-', 'open', false);
    st.mk = st.mk(iM);
    plot_ber_curve(gca, get_curve(S, tags13b{iM}), st, sprintf('M = %d', Mlist(iM)));
end
plot_ber_curve(gca, get_curve(S, 'f13b_ideal'), ...
    struct('col', [0 0 0], 'mk', 'none', 'ls', ':', 'open', true), 'точное положение кадра');
set(gca, 'YScale', 'log', 'XScale', 'linear');
xlim(gca, [cfgShort.snrLow(1)-0.5, cfgShort.snrLow(end)+0.5]);
ylim(gca, [0.1 0.6]); yticks(gca, [0.1 0.15 0.2 0.3 0.4 0.5]);
xticks(gca, cfgShort.snrLow(1):2:cfgShort.snrLow(end));
xlabel(gca, 'E_b/N_0, дБ'); ylabel(gca, 'BER');
title(gca, '(б) N_{CP} = 32: накопление помогает');
legend(gca, 'Location', 'southwest');
style_axes(gca);
save_figure(cfg, h, 'fig_1_3_ber_M');

% 1.4 Автокорреляционная функция преамбулы Минна и её метрика
h = new_figure(cfg, 'АКФ преамбулы Минна', 'Автокорреляционная функция преамбулы Минна');
tiledlayout(h, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
[acf, lags] = cyclic_acf(cfg.preambleTime);
plot(lags, acf, 'Color', [0 0.45 0.74], 'LineWidth', 1.6);
hold on; grid on;
idxPeak = findpeaks_simple(acf, lags);
plot(lags(idxPeak), acf(idxPeak), 'o', 'Color', [0.85 0.33 0.10], ...
     'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 6);
xlabel('сдвиг \tau, отсчёты'); ylabel('|R(\tau)| / R(0)');
ylim([-0.2 1.05]);
title('(а) циклическая АКФ преамбулы [A A -A -A]');
legend({'АКФ', 'максимумы'}, 'Location', 'northeast'); style_axes();

nexttile;
[mMinn, axMin] = minn_metric_shape(cfg, make_scen(), cfg.snrDemo, cfg.nAvgMetric);
plot(axMin, mMinn, 'Color', [0.49 0.18 0.56], 'LineWidth', 1.6);
grid on;
xlabel('отсчёты относительно начала тела преамбулы');
ylabel('M_{Minn}(d), отн. ед.');
title(sprintf('(б) метрика Минна, Eb/N0 = %d дБ', cfg.snrDemo));
style_axes();
save_figure(cfg, h, 'fig_1_4_minn_acf');

% 1.5 Помехоустойчивость при синхронизации по Минну
tags15 = cell(1, 2);
for iEq = 1:2
    scen = scen_set(make_scen(), 'deltaMode', 'uniform', 'timeSync', 'minn', ...
                    'useEq', (iEq == 1));
    tag = sprintf('f15_eq%d', iEq == 1);
    out = ber_sweep(cfg, scen, cfg.snrDb, sprintf('Минн, %s', ...
          ternary(iEq == 1, 'с экв.', 'без экв.')));
    S = add_curve(S, tag, out); tags15{iEq} = tag;
end

h = new_figure(cfg, 'BER при синхронизации по Минну', 'Синхронизация по Минну');
ax = axes(h); hold(ax, 'on'); grid(ax, 'on');
st1 = struct('col', [0 0.45 0.74], 'mk', 'o', 'ls', '-',  'open', false);
st2 = struct('col', [0.85 0.33 0.10], 'mk', 's', 'ls', '--', 'open', true);
plot_ber_curve(ax, get_curve(S, tags15{1}), st1, 'Минн, с эквалайзером');
plot_ber_curve(ax, get_curve(S, tags15{2}), st2, 'Минн, без эквалайзера');
plot_theory(ax, cfg);
finish_ber_axes(ax, cfg, 'Временная синхронизация по преамбуле Минна');
save_figure(cfg, h, 'fig_1_5_minn_ber');

%% ЧАСТЬ 2. Частотная синхронизация
fprintf('\n--- Часть 2. Частотная синхронизация -------------------------\n');

% 2.1 Дробная расстройка без частотной синхронизации
epsList = [0.1 0.3 0.5];
tags21 = cell(1, 2*numel(epsList));
for iE = 1:numel(epsList)
    for iEq = 1:2
        scen = scen_set(make_scen(), 'eps', epsList(iE), 'deltaMode', 'uniform', ...
                        'timeSync', 'genie', 'useEq', (iEq == 1));
        tag = sprintf('f21_e%d_eq%d', iE, iEq == 1);
        out = ber_sweep(cfg, scen, cfg.snrDb, sprintf('\\epsilon = %.1f, %s', ...
              epsList(iE), ternary(iEq == 1, 'с экв.', 'без экв.')));
        S = add_curve(S, tag, out); tags21{(iE-1)*2 + iEq} = tag;
    end
end

h = new_figure(cfg, 'BER при расстройке по частоте', 'Влияние дробной частотной расстройки');
ax = axes(h); hold(ax, 'on'); grid(ax, 'on');
colsE = [0.30 0.75 0.93; 0.93 0.69 0.13; 0.64 0.08 0.18];
for iE = 1:numel(epsList)
    for iEq = 1:2
        st = struct('col', colsE(iE,:), 'mk', 'os^', 'ls', ternary(iEq == 1, '-', '--'), ...
                    'open', iEq == 2);
        st.mk = st.mk(iE);
        plot_ber_curve(ax, get_curve(S, tags21{(iE-1)*2 + iEq}), st, ...
            sprintf('\\epsilon = %.1f, %s', epsList(iE), ...
            ternary(iEq == 1, 'с эквалайзером', 'без эквалайзера')));
    end
end
plot_theory(ax, cfg);
finish_ber_axes(ax, cfg, 'Без частотной синхронизации');
save_figure(cfg, h, 'fig_2_1_cfo');

% 2.2 Влияние M на частотную синхронизацию по ЦП
tags22 = cell(1, numel(Mlist));
for iM = 1:numel(Mlist)
    scen = scen_set(make_scen(), 'epsMode', 'uniform', 'epsRange', [-0.5 0.5], ...
                    'deltaMode', 'uniform', 'timeSync', 'cp', 'freqSync', 'cp', ...
                    'M', Mlist(iM), 'useEq', true);
    tag = sprintf('f22_M%d', Mlist(iM));
    out = ber_sweep(cfg, scen, cfg.snrDb, sprintf('ЦП-ЧС, M = %d', Mlist(iM)));
    S = add_curve(S, tag, out); tags22{iM} = tag;
end

h = new_figure(cfg, 'BER частотной синхронизации по ЦП', 'Частотная синхронизация по ЦП');
tiledlayout(h, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hold on; grid on;
for iM = 1:numel(Mlist)
    st = struct('col', cols(iM,:), 'mk', 'os^d', 'ls', '-', 'open', false);
    st.mk = st.mk(iM);
    plot_ber_curve(gca, get_curve(S, tags22{iM}), st, sprintf('M = %d', Mlist(iM)));
end
cfgRef = cfg;  cfgRef.nFrames = 40; cfgRef.minFrames = 30; cfgRef.maxErrors = Inf;
outNo22 = ber_sweep(cfgRef, scen_set(make_scen(), 'epsMode', 'uniform', ...
    'epsRange', [-0.5 0.5], 'deltaMode', 'uniform', 'timeSync', 'genie', ...
    'freqSync', 'none', 'useEq', true), cfg.snrDb, 'без синхр. частоты');
S = add_curve(S, 'f22_none', outNo22);
plot_ber_curve(gca, get_curve(S, 'f22_none'), ...
    struct('col', [0.5 0.5 0.5], 'mk', 'x', 'ls', ':', 'open', true), 'без синхронизации частоты');
plot_theory(gca, cfg);
finish_ber_axes(gca, cfg, '(а) BER, \epsilon \in (-0.5, 0.5)');

nexttile;
hold on; grid on;
CFOS = cell(1, numel(Mlist));
for iM = 1:numel(Mlist)
    CFOS{iM} = cfo_error_curve(cfg, [-0.5 0.5], Mlist(iM), cfg.snrErr, cfg.nFramesErr);
    fprintf('  СКО оценки eps (M=%d): %s\n', Mlist(iM), mat2str(CFOS{iM}, 2));
    mkList = 'os^d';
    loglog(gca, cfg.snrErr, CFOS{iM}, 'LineStyle', '-', 'Color', cols(iM,:), ...
           'Marker', mkList(iM), 'MarkerSize', 5, 'LineWidth', 1.4, ...
           'DisplayName', sprintf('M = %d', Mlist(iM)));
    % теория: sigma_eps = 1/(2*pi*sqrt(M*Ncp*gamma))
    gam = 10.^((cfg.snrErr + 10*log10(cfg.qamBits*numel(cfg.dataIdx)/cfg.Nused))/10);
    th = 1./(2*pi*sqrt(Mlist(iM)*cfg.Ncp*gam));
    loglog(gca, cfg.snrErr, th, '--', 'Color', cols(iM,:), 'LineWidth', 0.8, ...
           'HandleVisibility', 'off');
end
set(gca, 'XScale', 'linear', 'YScale', 'log');
xlabel(gca, 'E_b/N_0, дБ');
ylabel(gca, 'СКО оценки \epsilon');
title(gca, '(б) точность оценки ЧС (штрих - теория)');
legend(gca, 'Location', 'southwest');
style_axes(gca);
save_figure(cfg, h, 'fig_2_2_cp_freq_M');

% 2.3 Влияние длины ЦП
NcpList = [64 128 256 512];
tags23 = cell(1, numel(NcpList));
for iL = 1:numel(NcpList)
    cfgL = cfg;
    cfgL.Ncp  = NcpList(iL);
    cfgL.Nsym = cfgL.Nfft + cfgL.Ncp;
    scen = scen_set(make_scen(), 'epsMode', 'uniform', 'epsRange', [-0.5 0.5], ...
                    'deltaMode', 'uniform', 'timeSync', 'cp', 'freqSync', 'cp', ...
                    'M', 1, 'useEq', true);
    tag = sprintf('f23_Ncp%d', NcpList(iL));
    out = ber_sweep(cfgL, scen, cfg.snrDb, sprintf('N_{CP} = %d', NcpList(iL)));
    S = add_curve(S, tag, out); tags23{iL} = tag;
end

h = new_figure(cfg, 'BER против длины ЦП', 'Влияние длины циклического префикса');
ax = axes(h); hold(ax, 'on'); grid(ax, 'on');
colsL = [0.47 0.67 0.19; 0.93 0.69 0.13; 0.85 0.33 0.10; 0.49 0.18 0.56];
for iL = 1:numel(NcpList)
    st = struct('col', colsL(iL,:), 'mk', 'os^d', 'ls', '-', 'open', false);
    st.mk = st.mk(iL);
    plot_ber_curve(ax, get_curve(S, tags23{iL}), st, ...
        sprintf('N_{CP} = %d (%.0f%%)', NcpList(iL), 100*NcpList(iL)/cfg.Nfft));
end
plot_theory(ax, cfg);
finish_ber_axes(ax, cfg, 'Частотная синхронизация по ЦП, \epsilon \in (-0.5, 0.5)');
save_figure(cfg, h, 'fig_2_3_cp_len');

% 2.4 Частотная синхронизация по преамбуле Минна
epsRanges = {[-0.5 0.5], [-2 2]};
tags24 = cell(1, numel(epsRanges));
for iR = 1:numel(epsRanges)
    scen = scen_set(make_scen(), 'epsMode', 'uniform', 'epsRange', epsRanges{iR}, ...
                    'deltaMode', 'uniform', 'timeSync', 'minn', 'freqSync', 'minn', ...
                    'useEq', true);
    tag = sprintf('f24_r%d', iR);
    out = ber_sweep(cfg, scen, cfg.snrDb, sprintf('Минн-ЧС, \\epsilon \\in (%.1f, %.1f)', ...
          epsRanges{iR}(1), epsRanges{iR}(2)));
    S = add_curve(S, tag, out); tags24{iR} = tag;
end
% опорная кривая: расстройка без синхронизации
scenNo = scen_set(make_scen(), 'epsMode', 'uniform', 'epsRange', [-2 2], ...
                  'deltaMode', 'uniform', 'timeSync', 'genie', 'freqSync', 'none', ...
                  'useEq', true);
outRef = ber_sweep(cfgRef, scenNo, cfg.snrDb, 'без ЧС, \epsilon \in (-2, 2)');
S = add_curve(S, 'f24_ref', outRef);

h = new_figure(cfg, 'BER частотной синхронизации по Минну', 'Частотная синхронизация по преамбуле Минна');
ax = axes(h); hold(ax, 'on'); grid(ax, 'on');
st = struct('col', [0 0.45 0.74], 'mk', 'o', 'ls', '-', 'open', false);
plot_ber_curve(ax, get_curve(S, tags24{1}), st, '\epsilon \in (-0.5, 0.5)');
st = struct('col', [0.85 0.33 0.10], 'mk', 's', 'ls', '--', 'open', true);
plot_ber_curve(ax, get_curve(S, tags24{2}), st, '\epsilon \in (-2, 2)');
st = struct('col', [0.5 0.5 0.5], 'mk', 'x', 'ls', ':', 'open', true);
plot_ber_curve(ax, get_curve(S, 'f24_ref'), st, 'без синхронизации, \epsilon \in (-2, 2)');
plot_theory(ax, cfg);
finish_ber_axes(ax, cfg, 'Оценка по преамбуле [A A -A -A]');
save_figure(cfg, h, 'fig_2_4_minn_freq');


%% Сводка результатов
fprintf('\n=============================================================\n');
fprintf(' СВОДКА: BER при Eb/N0 = %g дБ\n', cfg.snrRef);
fprintf('-------------------------------------------------------------\n');
for i = 1:numel(S)
    [~, k] = min(abs(S(i).snr - cfg.snrRef));
    fprintf(' %-34s BER = %.3e  (%d ошибок / %d бит)\n', ...
        S(i).tag, S(i).ber(k), S(i).err(k), S(i).bits(k));
end
fprintf('-------------------------------------------------------------\n');
fprintf(' Полное время моделирования: %.1f с\n', toc(tAll));
fprintf(' Графики сохранены в: %s\n', cfg.figDir);

save(fullfile(cfg.outDir, 'Lab3_results.mat'), 'S', 'cfg');
diary off;


%% ЛОКАЛЬНЫЕ ФУНКЦИИ

function cfg = lab3_cfg()
%LAB3_CFG  Единая точка настройки моделирования.
    cfg.Nfft      = 2048;                 % размер БПФ (число поднесущих)
    cfg.Ncp       = 256;                  % длина ЦП (1/8 от Nfft)
    cfg.pilotStep = 4;                    % пилот на каждой 4-й поднесущей
    cfg.nDataSym  = 14;                   % OFDM-символов данных в кадре
    cfg.evalSym   = 3:12;                 % из них используются для оценки BER
    cfg.qamBits   = 4;                    % КАМ-16 -> 4 бита на символ
    cfg.padLen    = cfg.Nfft;             % запас нулей (для отрицательных сдвигов)
    cfg.snrDb     = 0:2:30;               % Eb/N0, дБ
    cfg.snrRef    = 14;                   % опорное Eb/N0 для сводной таблицы
    cfg.snrLow    = -12:2:0;              % низко-ОСШ диапазон (короткий ЦП)
    cfg.snrDemo   = 0;                    % Eb/N0 для демонстрационных рисунков
    cfg.nFrames   = 120;                  % максимум кадров на точку BER
    cfg.minFrames = 10;                   % минимум кадров (устойчивость оценки)
    cfg.maxErrors = 400;                  % ранняя остановка по числу ошибок
    cfg.nAvgMetric = 20;                  % усреднение метрики по кадрам
    cfg.snrDelta  = 12;                   % Eb/N0 для развёртки по delta
    cfg.nFramesDelta = 40;                % кадров на точку развёртки по delta
    cfg.snrErr    = -10:2:10;             % Eb/N0 для оценки точности ЧС
    cfg.nFramesErr = 20;                  % кадров на точку кривой точности
    cfg.seed      = 20241003;
    cfg.minChannelGain = 0.1;             % защита при делении на оценку канала
    cfg.figureResolution = 170;

    % быстрый отладочный режим
    if strcmpi(getenv('LAB3_QUICK'), '1')
        cfg.nFrames = 6; cfg.maxErrors = 40; cfg.snrDb = 0:6:12; cfg.nAvgMetric = 4;
        cfg.nFramesDelta = 4; cfg.snrLow = -12:6:0;
    end

    cfg.Nsym  = cfg.Nfft + cfg.Ncp;
    cfg.Nused = cfg.Nfft - 1;             % активные поднесущие (без DC)

    % --- карта поднесущих ---------------------------------------------
    kk = 0:cfg.Nfft-1;
    cfg.dcIdx    = cfg.Nfft/2 + 1;                          % нулевая (DC) поднесущая
    cfg.pilotIdx = find(mod(kk, cfg.pilotStep) == 0).';      % пилоты (КАМ-4)
    cfg.pilotIdx = setdiff(cfg.pilotIdx, cfg.dcIdx).';
    cfg.dataIdx  = setdiff((1:cfg.Nfft).', union(cfg.pilotIdx, cfg.dcIdx));

    % --- известные последовательности (детерминированные) -------------
    cfg.pilotSyms    = qpsk_seq(numel(cfg.pilotIdx), 12345);
    cfg.preambleTime = minn_preamble_time(cfg);              % [A A -A -A]
    cfg.preambleX    = fft(cfg.preambleTime) / sqrt(cfg.Nfft);   % спектр преамбулы

    % --- каталоги и видимость графиков --------------------------------
    p = fileparts(mfilename('fullpath'));
    if isempty(p); p = pwd; end
    cfg.outDir = p;
    cfg.figDir = fullfile(p, 'figures');
    if usejava('desktop'); cfg.figVisible = 'on'; else; cfg.figVisible = 'off'; end
end

% ------------------------------------------------------------------------
function [frame, tx] = lab_frame(cfg)
%LAB_FRAME  Формирование кадра: преамбула Минна + OFDM-символы с пилотами.
    N = cfg.Nfft;  Ncp = cfg.Ncp;
    nSym = 1 + cfg.nDataSym;
    X = zeros(N, nSym);
    X(:, 1) = cfg.preambleX;                       % преамбула (все поднесущие)

    % биты полезной нагрузки: 4 x Ndata x nDataSym
    bits = randi([0 1], cfg.qamBits, numel(cfg.dataIdx), cfg.nDataSym);
    syms = qam16_map(reshape(bits, cfg.qamBits, []));          % 1 x (Ndata*nSym)
    X(cfg.dataIdx, 2:end)  = reshape(syms, numel(cfg.dataIdx), cfg.nDataSym);
    X(cfg.pilotIdx, 2:end) = repmat(cfg.pilotSyms, 1, cfg.nDataSym);

    x = ifft(X, N, 1) * sqrt(N);                   % OFDM-модуляция (Parseval)
    xcp = [x(N-Ncp+1:N, :); x];                    % добавление ЦП
    frame = xcp(:);
    tx.bits = bits;
    tx.X = X;
end

% ------------------------------------------------------------------------
function [rx, ch] = lab_channel(cfg, frame, delta, epsNorm, ebnoDb)
%LAB_CHANNEL  Канал: временной сдвиг, частотный сдвиг, АБГШ.
%   Соглашение о знаке delta (как в методичке):
%     -Ncp <= delta < 0 - окно БПФ начинается внутри ЦП: МСИ нет, только
%                         фазовый наклон по поднесущим;
%     delta < -Ncp      - окно выходит за пределы ЦП: появляется МСИ;
%     delta > 0         - окно запаздывает: МСИ из следующего символа.
    N = cfg.Nfft;  P = cfg.padLen;
    nTail = P + abs(delta) + 8;
    rx = zeros(numel(frame) + P + nTail, 1);
    start = P + 1 - delta;                          % ПЕРВЫЙ отсчёт кадра
    rx(start:start+numel(frame)-1) = frame;

    % частотный сдвиг: r[n] = s[n] * exp(j*2*pi*eps*n/N)
    n = (0:numel(rx)-1).';
    if epsNorm ~= 0
        rx = rx .* exp(1j*2*pi*epsNorm*n/N);
    end

    % АБГШ: мощность шума из Eb/N0 с учётом служебных поднесущих
    gammaDb = ebnoDb + 10*log10(cfg.qamBits*numel(cfg.dataIdx)/cfg.Nused);
    sigma2  = (cfg.Nused/N) / (10^(gammaDb/10));
    rx = rx + sqrt(sigma2/2)*(randn(size(rx)) + 1j*randn(size(rx)));

    ch.delta = delta; ch.eps = epsNorm; ch.sigma2 = sigma2;
end

% ------------------------------------------------------------------------
function [dHat, info] = lab_sync(cfg, scen, rx, delta)
%LAB_SYNC  Временная и частотная синхронизация приёмника.
%   dHat - оценка положения ПЕРВОГО отсчёта ЦП первого символа кадра.
    dNom = cfg.padLen + 1;                          % ожидаемое начало кадра
    info.metric = []; info.axis = [];

    switch scen.timeSync
        case 'none'                                 % синхронизация не выполняется
            dHat = dNom;
        case 'genie'                                % точное положение (идеал)
            dHat = dNom - delta;
        case 'cp'                                   % по циклическому префиксу
            [dHat, info.metric, info.axis] = cp_metric(rx, cfg, scen.M, dNom);
        case 'minn'                                 % по преамбуле Минна
            [dHat, info.metric, info.axis] = minn_metric(rx, cfg, dNom);
        otherwise
            error('lab_sync:badMode', 'Неизвестный метод: %s', scen.timeSync);
    end

    switch scen.freqSync
        case 'none'
            info.epsHat = 0;
        case 'cp'
            info.epsHat = cfo_from_cp(rx, cfg, dHat, scen.M);
        case 'minn'
            info.epsHat = cfo_from_minn(rx, cfg, dHat);
        otherwise
            error('lab_sync:badMode', 'Неизвестный метод ЧС: %s', scen.freqSync);
    end
end

% ------------------------------------------------------------------------
function [nErr, nBit] = frame_trial(cfg, scen, ebnoDb)
%FRAME_TRIAL  Один кадр: передача, канал, синхронизация, демодуляция, BER.
    [txFrame, tx] = lab_frame(cfg);

    [delta, epsNorm] = scen_impairment(cfg, scen);
    rx = lab_channel(cfg, txFrame, delta, epsNorm, ebnoDb);

    [dHat, info] = lab_sync(cfg, scen, rx, delta);
    rxc = cfo_correct(cfg, rx, info.epsHat);

    % окна БПФ для выбранных символов кадра
    kSym = 1 + cfg.evalSym;                         % сдвиг на преамбулу
    win  = dHat + cfg.Ncp + (kSym - 1)*cfg.Nsym;
    idx  = win + (0:cfg.Nfft-1).';
    if any(idx(:) < 1) || any(idx(:) > numel(rxc))
        nErr = 0; nBit = 0; return;
    end

    Y = fft(rxc(idx), [], 1) / sqrt(cfg.Nfft);      % N x nEval

    if scen.useEq
        H = lab_channel_estimate(cfg, Y);
        Y = Y ./ H;
    end

    bitsHat = qam16_demap(Y(cfg.dataIdx, :));
    bitsRef = tx.bits(:, :, cfg.evalSym);
    nErr = sum(bitsHat(:) ~= bitsRef(:));
    nBit = numel(bitsRef);
end

% ------------------------------------------------------------------------
function [delta, epsNorm] = scen_impairment(cfg, scen)
%SCEN_IMPAIRMENT  Параметры расстройки для текущего кадра.
    switch scen.deltaMode
        case 'fixed'
            delta = scen.delta;
        case 'uniform'
            delta = randi([-cfg.Nfft/2, cfg.Nfft/2 - 1]);
        otherwise
            error('scen:badDelta', 'Неизвестный режим delta');
    end
    switch scen.epsMode
        case 'fixed'
            epsNorm = scen.eps;
        case 'uniform'
            epsNorm = scen.epsRange(1) + diff(scen.epsRange)*rand;
        otherwise
            error('scen:badEps', 'Неизвестный режим eps');
    end
end

% ------------------------------------------------------------------------
function H = lab_channel_estimate(cfg, Y)
%LAB_CHANNEL_ESTIMATE  Оценка канала по пилотам (два прохода + наклон ФЧХ).
%   1) МНК-оценка на пилотных поднесущих: Hp = Y(p)/P(p);
%   2) "форма" Hn = Hp/Hp(1,:) одинакова для всех символов кадра (в ней
%      сохраняется фазовый наклон от временного сдвига), поэтому усредняется
%      по символам - шум оценки падает в sqrt(nSym) раз;
%   3) по среднему приращению фазы между пилотами оценивается наклон ФЧХ
%      (остаточная задержка) и форма "выпрямляется" перед интерполяцией -
%      иначе линейная интерполяция комплексных чисел на крутом наклоне даёт
%      большую ошибку амплитуды;
%   4) общая фаза каждого символа cpe(m) оценивается как среднее Hp/Hshape по
%      всем 511 пилотам (шум падает в sqrt(511) раз);
%   5) H(k,m) = Hshape(k) * cpe(m).
    pos  = cfg.pilotIdx(:);
    grid = (1:cfg.Nfft).';
    Hp   = Y(cfg.pilotIdx, :) ./ cfg.pilotSyms;     % 511 x nSym, МНК на пилотах

    ref = Hp(1, :);                                % опорный пилот каждого символа
    ref(abs(ref) < realmin) = realmin;
    Hn  = Hp ./ ref;                               % форма, общая для всех символов
    Hm  = mean(Hn, 2);                             % 511 x 1, усреднение по символам

    dph   = angle(Hm(2:end) .* conj(Hm(1:end-1))); % приращения фазы между пилотами
    slope = median(dph ./ diff(pos));              % рад на поднесущую (наклон ФЧХ)
    Hd    = Hm .* exp(-1j*slope*pos);              % "выпрямленная" форма

    Hr = interp1(pos, real(Hd), grid, 'linear', 'extrap');
    Hi = interp1(pos, imag(Hd), grid, 'linear', 'extrap');
    Hshape = (Hr + 1j*Hi) .* exp(1j*slope*grid);   % 2048 x 1

    cpe = mean(Hp ./ Hm, 1);                       % 1 x nSym, общая фаза символа
    H   = Hshape * cpe;                            % 2048 x nSym

    mag = abs(H);
    bad = mag < cfg.minChannelGain;                % защита от деления на ~0
    H(bad) = cfg.minChannelGain * exp(1j*angle(H(bad)));
end

% ------------------------------------------------------------------------
function [dHat, metric, axisD] = cp_metric(rx, cfg, M, dNom)
%CP_METRIC  Метрика временной синхронизации по ЦП с когерентным накоплением.
%   M(d) = |sum_m A(d + m*Nsym)|^2 / (sum_m P(d + m*Nsym))^2,
%   A(d) = sum_{i=0}^{Ncp-1} r*(d+i) r(d+i+Nfft)  - корреляция ЦП и его копии.
    N = cfg.Nfft;  L = cfg.Ncp;  Nsym = cfg.Nsym;
    r = rx(:);
    r2 = abs(r).^2;

    c  = conj(r(1:end-N)) .* r(1+N:end);            % попарные произведения
    cs = [0; cumsum(c)];                            % префиксные суммы
    e  = [0; cumsum(r2)];

    nc   = numel(c);
    dmax = nc - L + 1;
    dd   = (1:dmax).';
    A    = cs(dd+L) - cs(dd);                                       % числитель
    Pn   = (e(dd+L) - e(dd)) + (e(dd+N+L) - e(dd+N));               % энергия

    dValid = dmax - (M-1)*Nsym;
    if dValid < 1; error('cp_metric:short', 'Сигнал короче окна накопления'); end

    Aacc = A(1:dValid);
    Pacc = Pn(1:dValid);
    for m = 1:M-1
        s = m*Nsym;
        Aacc = Aacc + A(1+s : dValid+s);
        Pacc = Pacc + Pn(1+s : dValid+s);
    end
    metric = abs(Aacc).^2 ./ max(Pacc.^2, realmin);

    lo = max(1, dNom - floor(N/2) - 64);
    hi = min(dValid, dNom + ceil(N/2) + 64);
    [~, k] = max(metric(lo:hi));
    dHat   = lo + k - 1;                            % метрика ЦП даёт начало ЦП
    axisD  = (lo:hi).' - dNom;
    metric = metric(lo:hi);
end

% ------------------------------------------------------------------------
function [dHat, metric, axisD] = minn_metric(rx, cfg, dNom)
%MINN_METRIC  Метрика Минна по преамбуле [A A -A -A].
%   P(d) = sum_{m=0}^{2} b(m) * sum_{n} r*(d+n+mL) r(d+n+(m+1)L),
%   b = [+1 -1 +1]  - веса Минна, устраняющие "плато" метрики Шмидла-Кокса.
    N = cfg.Nfft;  L = N/4;
    r = rx(:);

    c  = conj(r(1:end-L)) .* r(1+L:end);
    cs = [0; cumsum(c)];
    e  = [0; cumsum(abs(r).^2)];

    nc   = numel(c);
    dmax = nc - L + 1;
    dd   = (1:dmax).';
    P0   = cs(dd+L) - cs(dd);                       % корреляция соседних блоков

    dValid = numel(r) - 4*L + 1;
    if dValid < 1; error('minn_metric:short', 'Сигнал короче преамбулы'); end
    P = P0(1:dValid) - P0(1+L : dValid+L) + P0(1+2*L : dValid+2*L);
    R = e((1:dValid) + 4*L) - e(1:dValid);          % энергия четырёх блоков
    metric = abs(P).^2 ./ max(R.^2, realmin);

    lo = max(1, dNom + cfg.Ncp - floor(N/2) - 64);
    hi = min(dValid, dNom + cfg.Ncp + ceil(N/2) + 64);
    [~, k] = max(metric(lo:hi));
    dHat   = max(1, lo + k - 1 - cfg.Ncp);          % метрика Минна даёт начало тела
    axisD  = (lo:hi).' - (dNom + cfg.Ncp);
    metric = metric(lo:hi);
end

% ------------------------------------------------------------------------
function epsHat = cfo_from_cp(rx, cfg, dHat, M)
%CFO_FROM_CP  Оценка дробной частотной расстройки по ЦП (диапазон |eps|<0.5).
    N = cfg.Nfft;  L = cfg.Ncp;  Nsym = cfg.Nsym;
    r = rx(:);
    acc = 0;
    for m = 0:M-1
        i0 = dHat + m*Nsym;
        acc = acc + sum(conj(r(i0 + (0:L-1))) .* r(i0 + N + (0:L-1)));
    end
    epsHat = angle(acc) / (2*pi);
end

% ------------------------------------------------------------------------
function epsHat = cfo_from_minn(rx, cfg, dHat)
%CFO_FROM_MINN  Оценка расстройки по преамбуле Минна (диапазон |eps|<2).
%   dHat - начало ЦП преамбулы, поэтому блоки начинаются с dHat + Ncp.
    N = cfg.Nfft;  L = N/4;
    r = rx(:);
    i0 = dHat + cfg.Ncp;                            % начало тела преамбулы
    p0 = sum(conj(r(i0 + (0:L-1)))       .* r(i0 + L   + (0:L-1)));
    p1 = sum(conj(r(i0 + L + (0:L-1)))   .* r(i0 + 2*L + (0:L-1)));
    p2 = sum(conj(r(i0 + 2*L + (0:L-1))) .* r(i0 + 3*L + (0:L-1)));
    acc = p0 - p1 + p2;                             % веса Минна [+1 -1 +1]
    epsHat = angle(acc) * N / (2*pi*L);
end

% ------------------------------------------------------------------------
function rxc = cfo_correct(cfg, rx, epsHat)
%CFO_CORRECT  Компенсация частотной расстройки.
    if epsHat == 0
        rxc = rx;
        return;
    end
    n = (0:numel(rx)-1).';
    rxc = rx .* exp(-1j*2*pi*epsHat*n/cfg.Nfft);
end

% ------------------------------------------------------------------------
function s = qam16_map(bits)
%QAM16_MAP  КАМ-16 с кодами Грея, средняя мощность 1.
%   bits - матрица 4 x N (строки: b1 b2 -> I, b3 b4 -> Q).
    b1 = bits(1,:); b2 = bits(2,:); b3 = bits(3,:); b4 = bits(4,:);
    vi = 2*b1 + xor(b1, b2);                        % 0..3 (00,01,11,10)
    vq = 2*b3 + xor(b3, b4);
    s  = ((2*vi - 3) + 1j*(2*vq - 3)) / sqrt(10);
end

% ------------------------------------------------------------------------
function bits = qam16_demap(syms)
%QAM16_DEMAP  Жёсткие решения по КАМ-16 (обратное к QAM16_MAP).
    lev = [-3 -1 1 3];
    v   = sqrt(10) * syms(:);                       % возврат к уровням +-1,+-3
    [~, vi] = min(abs(real(v) - lev), [], 2);
    [~, vq] = min(abs(imag(v) - lev), [], 2);
    vi = vi - 1;  vq = vq - 1;
    b1 = vi >= 2;  b2 = (vi == 1) | (vi == 2);
    b3 = vq >= 2;  b4 = (vq == 1) | (vq == 2);
    bits = [b1, b2, b3, b4].';                  % 4 x N, как у qam16_map
end

% ------------------------------------------------------------------------
function s = qpsk_seq(n, seed)
%QPSK_SEQ  Детерминированная последовательность КАМ-4 единичной мощности.
%   Генератор Лемера (Парка-Миллера): x <- 48271*x mod (2^31-1).
%   Все произведения < 2^53, поэтому арифметика в double точна.
    a = 48271;  m = 2147483647;
    x = mod(abs(seed), m);
    if x == 0; x = 1; end
    s = zeros(n, 1);
    for i = 1:n
        x = mod(a*x, m);
        s(i) = ((1 - 2*mod(x, 2)) + 1j*(1 - 2*mod(floor(x/2), 2))) / sqrt(2);
    end
end

% ------------------------------------------------------------------------
function y = minn_preamble_time(cfg)
%MINN_PREAMBLE_TIME  Преамбула [A A -A -A] (4 одинаковых блока, знаки + + - -).
    N = cfg.Nfft;  L = N/4;
    q = qpsk_seq(L, 54321);
    A = ifft(q, L) * sqrt(L);                       % блок единичной мощности
    y = [A; A; -A; -A];
    y = y / sqrt(mean(abs(y).^2));                  % нормировка мощности
end

% ------------------------------------------------------------------------
function [acf, lags] = cyclic_acf(y)
%CYCLIC_ACF  Циклическая автокорреляционная функция (через спектр).
    y = y(:);
    N = numel(y);
    Y = fft(y);
    R = ifft(abs(Y).^2) / N;
    lags = (-N/2:N/2-1).';
    acf  = abs(fftshift(R));
    acf  = acf / max(acf);
end

% ------------------------------------------------------------------------
function idx = findpeaks_simple(acf, lags)
%FINDPEAKS_SIMPLE  Локальные максимумы АКФ (без Signal Processing Toolbox).
    idx = [];
    for i = 2:numel(acf)-1
        if acf(i) > acf(i-1) && acf(i) >= acf(i+1) && acf(i) > 0.2
            idx(end+1) = i; %#ok<AGROW>
        end
    end
    idx = idx(:);
    [~, order] = sort(abs(lags(idx)));
    idx = idx(order(1:min(6, numel(order))));
end

% ------------------------------------------------------------------------
function [mAvg, axisD] = metric_shape(cfg, scen, M, ebnoDb, nAvg)
%METRIC_SHAPE  Усреднённая по кадрам метрика временной синхронизации по ЦП.
    mAvg = [];  axisD = [];
    for i = 1:nAvg
        txFrame = lab_frame(cfg);
        [delta, epsNorm] = scen_impairment(cfg, scen);
        rx = lab_channel(cfg, txFrame, delta, epsNorm, ebnoDb);
        [~, m, ax] = cp_metric(rx, cfg, M, cfg.padLen + 1);
        if isempty(mAvg)
            mAvg = zeros(size(m));  axisD = ax;
        end
        mAvg = mAvg + m;
    end
    mAvg = mAvg / nAvg;
end

% ------------------------------------------------------------------------
function [mAvg, axisD] = minn_metric_shape(cfg, scen, ebnoDb, nAvg)
%MINN_METRIC_SHAPE  Усреднённая по кадрам метрика Минна.
    mAvg = []; axisD = [];
    for i = 1:nAvg
        txFrame = lab_frame(cfg);
        [delta, epsNorm] = scen_impairment(cfg, scen);
        rx = lab_channel(cfg, txFrame, delta, epsNorm, ebnoDb);
        [~, m, ax] = minn_metric(rx, cfg, cfg.padLen + 1);
        if isempty(mAvg); mAvg = zeros(size(m)); axisD = ax; end
        mAvg = mAvg + m;
    end
    mAvg = mAvg / nAvg;
end

% ------------------------------------------------------------------------
function rmsErr = cfo_error_curve(cfg, epsRange, M, snrDb, nFrames)
%CFO_ERROR_CURVE  СКО оценки нормированной ЧС при точном положении кадра.
    rng(cfg.seed);
    rmsErr = zeros(size(snrDb));
    for i = 1:numel(snrDb)
        acc = 0;
        for f = 1:nFrames
            txFrame = lab_frame(cfg);
            epsTrue = epsRange(1) + diff(epsRange)*rand;
            delta   = randi([-cfg.Nfft/2, cfg.Nfft/2 - 1]);
            rx = lab_channel(cfg, txFrame, delta, epsTrue, snrDb(i));
            dHat = cfg.padLen + 1 - delta;          % точное положение кадра
            epsHat = cfo_from_cp(rx, cfg, dHat, M);
            % ошибка считается "по кругу": оценка определена с точностью до
            % целого числа периодов, иначе редкие переходы через +-0.5 дают
            % выбросы, которые портят СКО
            err = mod(epsHat - epsTrue + 0.5, 1) - 0.5;
            acc = acc + err^2;
        end
        rmsErr(i) = sqrt(acc/nFrames);
    end
end

% ------------------------------------------------------------------------
function ber = ber_point(cfg, scen, ebnoDb, nFrames, maxErr)
%BERR_POINT  BER при фиксированном Eb/N0 (для развёрток по параметру).
    rng(cfg.seed);
    nErr = 0; nBit = 0; nf = 0;
    while nf < nFrames && nErr < maxErr
        nf = nf + 1;
        [e, n] = frame_trial(cfg, scen, ebnoDb);
        nErr = nErr + e;  nBit = nBit + n;
    end
    if nBit > 0
        ber = nErr / nBit;
    else
        ber = NaN;
    end
    if ber == 0; ber = NaN; end
end

% ------------------------------------------------------------------------
function out = ber_sweep(cfg, scen, snrDb, label)
%BER_SWEEP  Монте-Карло: BER от Eb/N0 с ранней остановкой по ошибкам.
    rng(cfg.seed);                                  % общий поток случайных чисел
    out.snr = snrDb(:).';
    out.ber = nan(size(snrDb));
    out.err = zeros(size(snrDb));
    out.bits = zeros(size(snrDb));
    out.label = label;

    for i = 1:numel(snrDb)
        nErr = 0; nBit = 0; nFrames = 0;
        while nFrames < cfg.nFrames && (nErr < cfg.maxErrors || nFrames < cfg.minFrames)
            nFrames = nFrames + 1;
            [e, n] = frame_trial(cfg, scen, snrDb(i));
            nErr = nErr + e;  nBit = nBit + n;
        end
        out.err(i) = nErr;  out.bits(i) = nBit;
        if nBit > 0; out.ber(i) = nErr/nBit; end
        fprintf('  %-40s Eb/N0=%4.0f дБ: BER=%.3e (%4d/%6d, кадров %3d)\n', ...
                label, snrDb(i), out.ber(i), nErr, nBit, nFrames);
    end
end

% ------------------------------------------------------------------------
function s = make_scen()
%MAKE_SCEN  Сценарий по умолчанию.
    s = struct('deltaMode', 'fixed', 'delta', 0, ...
               'epsMode', 'fixed', 'eps', 0, 'epsRange', [-0.5 0.5], ...
               'timeSync', 'genie', 'freqSync', 'none', ...
               'M', 1, 'useEq', true);
end

% ------------------------------------------------------------------------
function s = scen_set(s, varargin)
%SCEN_SET  Изменение полей сценария парами 'имя', значение.
    for i = 1:2:numel(varargin)
        s.(varargin{i}) = varargin{i+1};
    end
end

% ------------------------------------------------------------------------
function S = add_curve(S, tag, out)
%ADD_CURVE  Добавление рассчитанной кривой в общий список.
    entry.tag  = tag;
    entry.snr  = out.snr;
    entry.ber  = out.ber;
    entry.err  = out.err;
    entry.bits = out.bits;
    if isempty(S); S = entry; else; S(end+1) = entry; end
end

% ------------------------------------------------------------------------
function c = get_curve(S, tag)
%GET_CURVE  Поиск кривой по метке.
    idx = find(strcmp({S.tag}, tag), 1);
    if isempty(idx)
        error('get_curve:notFound', 'Нет кривой с меткой %s', tag);
    end
    c = S(idx);
end

% ------------------------------------------------------------------------
function plot_ber_curve(ax, c, st, name)
%PLOT_BER_CURVE  Кривая BER (масштаб оси задаётся в finish_ber_axes).
    y = c.ber(:).';
    y(y <= 0) = NaN;
    pl = plot(ax, c.snr(:).', y, 'LineStyle', st.ls, 'Color', st.col, ...
        'Marker', st.mk, 'MarkerSize', 6, 'LineWidth', 1.4, 'DisplayName', name);
    if st.open; pl.MarkerFaceColor = 'none'; else; pl.MarkerFaceColor = st.col; end
end

% ------------------------------------------------------------------------
function plot_theory(ax, cfg)
%PLOT_THEORY  Теоретическая кривая КАМ-16 с кодами Грея в канале с АБГШ.
%   Точная вероятность ошибки на бит для 16-КАМ с кодом Грея:
%       Pb = 0.75*Q(a) + 0.25*Q(3a),   a = sqrt(gamma/5),
%   где gamma = Es/N0 на поднесущую. Первое слагаемое - ошибка "в соседний
%   уровень", второе - редкая ошибка "через уровень". При больших gamma
%   второе слагаемое исчезает, и Pb = 0.75*Q(a) = SER/4, где
%   SER = 3*Q(a) - 2.25*Q(a)^2 - вероятность ошибки на символ.
    snr = linspace(cfg.snrDb(1), cfg.snrDb(end), 301);
    gammaDb = snr + 10*log10(cfg.qamBits*numel(cfg.dataIdx)/cfg.Nused);
    gam = 10.^(gammaDb/10);
    q = @(x) 0.5*erfc(x/sqrt(2));
    a = sqrt(gam/5);
    ber = 0.75*q(a) + 0.25*q(3*a);
    plot(ax, snr, ber, 'k-.', 'LineWidth', 1.3, 'DisplayName', 'теория, КАМ-16 (АБГШ)');
end

% ------------------------------------------------------------------------
function finish_ber_axes(ax, cfg, ttl)
%FINISH_BER_AXES  Оформление осей BER (логарифмический масштаб по BER).
%   Масштаб задаётся ЯВНО: при включённом hold on функции semilogy/loglog
%   НЕ переключают масштаб уже существующих осей.
    set(ax, 'YScale', 'log', 'XScale', 'linear');
    xlim(ax, [cfg.snrDb(1)-0.5, cfg.snrDb(end)+0.5]);
    ylim(ax, [1e-6 1]);
    yticks(ax, 10.^(-6:0));
    xticks(ax, cfg.snrDb(1):2:cfg.snrDb(end));
    xlabel(ax, 'E_b/N_0, дБ');
    ylabel(ax, 'BER');
    title(ax, ttl);
    legend(ax, 'Location', 'southwest', 'NumColumns', 1);
    grid(ax, 'on');
    style_axes(ax);
end

% ------------------------------------------------------------------------
function style_axes(ax)
%STYLE_AXES  Единое оформление осей.
    if nargin < 1; ax = gca; end
    grid(ax, 'on'); box(ax, 'on');
    ax.GridAlpha = 0.22; ax.MinorGridAlpha = 0.12;
    ax.FontName = 'Arial'; ax.FontSize = 11;
    set(ax, 'YMinorGrid', 'on', 'XMinorGrid', 'on');
end

% ------------------------------------------------------------------------
function h = new_figure(cfg, name, ttl)
%NEW_FIGURE  Создание окна рисунка.
    h = figure('Name', name, 'NumberTitle', 'off', 'Color', 'w', ...
               'Visible', cfg.figVisible, 'Position', [80 60 1000 630]);
    if ~isempty(ttl); set(h, 'Name', ttl); end
end

% ------------------------------------------------------------------------
function save_figure(cfg, h, name)
%SAVE_FIGURE  Сохранение рисунка в PNG.
    f = fullfile(cfg.figDir, [name '.png']);
    exportgraphics(h, f, 'Resolution', cfg.figureResolution);
    fprintf('  -> сохранён %s.png\n', name);
end

% ------------------------------------------------------------------------
function out = ternary(cond, a, b)
%TERNARY  Короткая запись условного оператора.
    if cond; out = a; else; out = b; end
end

% ------------------------------------------------------------------------
function lab3_debug(cfg)
%LAB3_DEBUG  Самопроверка модели (вызывается при LAB3_DEBUG=1).
%   Проверяет: структуру ЦП, корреляцию ЦП и его копии, метрику при точном
%   положении, согласованность модулятора/демодулятора КАМ-16.
    fid = fopen(fullfile(cfg.outDir, 'Lab3_selfcheck.txt'), 'w', 'n', 'UTF-8');
    frame = lab_frame(cfg);
    d0 = cfg.padLen + 1;
    N = cfg.Nfft; Ncp = cfg.Ncp; Nsym = cfg.Nsym;
    say = @(varargin) fprintf(fid, varargin{:});

    say('длина кадра: %d отсчётов (%.2f символа)\n', numel(frame), numel(frame)/Nsym);

    % 1) ЦП - точная копия хвоста символа
    for s = 1:3
        off = (s-1)*Nsym;
        a = frame(off + (1:Ncp));
        b = frame(off + N + (1:Ncp));
        rho = abs(sum(conj(a).*b)) / sqrt(sum(abs(a).^2)*sum(abs(b).^2));
        say('символ %d: rho(ЦП, хвост) = %.6f, мощность = %.4f\n', ...
            s, rho, mean(abs(frame(off + (1:Nsym))).^2));
    end

    % 2) модулятор/демодулятор КАМ-16
    bt = double(randi([0 1], 4, 1000));
    bh = qam16_demap(qam16_map(bt));
    say('КАМ-16: BER(map->demap) = %g (должно быть 0)\n', mean(bt(:) ~= bh(:)));

    % 3) метрика по ЦП при точном положении и без шума
    rx = lab_channel(cfg, frame, 0, 0, 100);
    [dHat, metric, ax] = cp_metric(rx, cfg, 1, d0);
    say('метрика ЦП: dHat = %d (истина %d), M(0) = %.4f (теория 0.25)\n', ...
        dHat, d0, metric(ax == 0));
    [dHatM, metM] = minn_metric(rx, cfg, d0);
    say('метрика Минна: dHat = %d (истина %d), max M = %.4f\n', ...
        dHatM, d0, max(metM));

    % 4) нормировка шума
    gDb = 10 + 10*log10(cfg.qamBits*numel(cfg.dataIdx)/cfg.Nused);
    say('Eb/N0 = 10 дБ -> Es/N0 = %.2f дБ, sigma^2 = %.5f\n', ...
        gDb, (cfg.Nused/N)/(10^(gDb/10)));
    fclose(fid);
    fprintf('Самопроверка записана в %s\n', fullfile(cfg.outDir, 'Lab3_selfcheck.txt'));
end




% ------------------------------------------------------------------------
function rxShifted = apply_fractional_delay(rx, deltaFrac)
%APPLY_FRACTIONAL_DELAY  Дробная задержка сигнала через ДПФ.
    if deltaFrac == 0
        rxShifted = rx;
        return;
    end
    N = numel(rx);
    k = [0:ceil(N/2)-1, -floor(N/2):-1].';
    f = k / N;
    phaseShift = exp(-1j * 2 * pi * f * deltaFrac);
    rxShifted = ifft(fft(rx) .* phaseShift);
end

% ------------------------------------------------------------------------
function [hChannel, Hfreq] = generate_multipath_channel(Nfft, maxDelay, numTaps)
%GENERATE_MULTIPATH_CHANNEL  Канал Рэлея с экспоненциальным профилем задержек.
    if maxDelay < 1; maxDelay = 1; end
    delays = round(linspace(0, maxDelay, numTaps));
    pdp = exp(-delays / max(maxDelay / 3, 1));
    pdp = pdp / sum(pdp);
    taps = sqrt(pdp/2) .* (randn(1, numTaps) + 1j*randn(1, numTaps));
    hChannel = zeros(maxDelay + 1, 1);
    hChannel(delays + 1) = taps;
    hChannel = hChannel / norm(hChannel);
    Hfreq = fft(hChannel, Nfft);
end


