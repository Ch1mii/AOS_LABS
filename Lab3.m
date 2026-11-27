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

fprintf(' ЛАБОРАТОРНАЯ РАБОТА №3: синхронизация OFDM в канале с АБГШ\n');
fprintf(' Nfft = %d, Ncp = %d (%.0f%% от Nfft), защитные полосы: %d нулей с каждого края\n', ...
    cfg.Nfft, cfg.Ncp, 100*cfg.Ncp/cfg.Nfft, cfg.Nguard);
fprintf(' Активных поднесущих: %d (данные: КАМ-16 на %d, пилоты: КАМ-4 на %d, DC: 1)\n', ...
    cfg.Nused, numel(cfg.dataIdx), numel(cfg.pilotIdx));
fprintf(' Кадр: преамбула Минна + %d OFDM-символов (оценка BER по символам %d...%d)\n', ...
    cfg.nDataSym, cfg.evalSym(1), cfg.evalSym(end));
fprintf(' Информационных бит на кадр: %d\n', ...
    cfg.qamBits * numel(cfg.dataIdx) * numel(cfg.evalSym));
fprintf(' Eb/N0 = %s дБ, фиксированное число кадров на точку: %d\n\n', ...
    mat2str(cfg.snrDb), cfg.nFrames);

S = struct('tag', {}, 'snr', {}, 'ber', {}, 'err', {}, 'bits', {});

cols = [0 0.45 0.74; 0.85 0.33 0.10; 0.49 0.18 0.56; 0.47 0.67 0.19];

%%  ЧАСТЬ 0. Демонстрация временного и частотного сдвига

fprintf('--- Часть 0. Временной и частотный сдвиг ---------------------\n');

rng(cfg.seed);
txFrame = lab_frame(cfg);

% (а) временной сдвиг: кадр запаздывает на delta отсчётов
deltaDemo = 300;
rxRef = lab_channel(cfg, txFrame, 0,          0, cfg.snrDemo);
rxDel = lab_channel(cfg, txFrame, deltaDemo,  0, cfg.snrDemo);

% (б,в) фазовая ошибка от временного сдвига внутри ЦП: до и после эквалайзера
%       созвездие показываем при высоком ОСШ, чтобы фазовый поворот был нагляден
deltaCP  = -cfg.Ncp/8;
snrConst = 22;                                  % Eb/N0 для демонстрации созвездий
rxCP = lab_channel(cfg, txFrame, deltaCP, 0, snrConst);
win0 = cfg.padLen + 1 + cfg.Ncp + (cfg.evalSym(1))*cfg.Nsym;
Seg  = rxCP(win0 + (0:cfg.Nfft-1).');
Ycp  = fftshift(fft(Seg));
Hcp  = lab_channel_estimate(cfg, Ycp);
YcpEq = Ycp ./ Hcp;

h = new_figure(cfg, 'Сдвиги сигнала', 'Временной и частотный сдвиг OFDM-сигнала');
tiledlayout(h, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% (а) временной сдвиг виден как смещение максимума корреляции по ЦП
[~, mRef, axMet] = cp_metric(rxRef, cfg, 1, cfg.padLen + 1, 'rho');
[~, mDel]        = cp_metric(rxDel, cfg, 1, cfg.padLen + 1, 'rho');
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
%     delta = 0 | -Ncp/8 (внутри ЦП) | -(Ncp+64) (за пределами ЦП) | +96
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
        berSweep(iD, iEq) = ber_point(cfg, scen, cfg.snrDelta, cfg.nFramesDelta);
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

% 1.2 Влияние числа символов когерентного накопления M на метрику (без нормировки на max)
Mlist = [1 2 4 8];
h = new_figure(cfg, 'Метрика ЦП', 'Автокорреляционная функция по ЦП');
ax = axes(h); hold(ax, 'on'); grid(ax, 'on');
mStyles = {'-', '--', '-.', ':'};
hleg = gobjects(numel(Mlist), 1);
maxPeak12 = 0;
for iM = 1:numel(Mlist)
    [mAvg, axD] = metric_shape(cfg, make_scen(), Mlist(iM), cfg.snrMetric12, cfg.nAvgMetric12);
    maxPeak12 = max(maxPeak12, max(mAvg));
    pl = plot(ax, axD, mAvg, 'LineStyle', mStyles{iM}, 'Color', cols(iM,:), 'LineWidth', 1.5, ...
              'DisplayName', sprintf('M = %d', Mlist(iM)));
    hleg(iM) = pl;
    far = abs(axD) > 1.1*cfg.Ncp;
    fprintf(['  M = %d: пик |АКФ| = %.5f, средний шумовой уровень вне пика ' ...
             '%.6f\n'], Mlist(iM), max(mAvg), mean(mAvg(far), 'omitnan'));
end
set(ax, 'YScale', 'linear');
xlim(ax, [-1.3*cfg.Ncp, 1.3*cfg.Ncp]);
ylim(ax, [0, maxPeak12 * 1.15]);
xlabel(ax, 'отсчёты задержки относительно истинного начала символа');
ylabel(ax, '|АКФ|');
title(ax, sprintf('АКФ по ЦП без нормировки на максимум (E_b/N_0 = %d дБ, %d реализаций)', ...
    cfg.snrMetric12, cfg.nAvgMetric12));
legend(hleg, 'Location', 'northeast'); style_axes(ax);
save_figure(cfg, h, 'fig_1_2_metric_M');

% 1.3 Влияние M на помехоустойчивость (слепая синхронизация по ЦП без подглядывания delta)
%     (а) БЕЗ эквалайзера (штриховые линии): ошибка синхронизации по ЦП вызывает
%         фазовый наклон exp(-j*2*pi*k*delta/Nfft), разрушающий КАМ-16.
%         Увеличение M снижает вероятность промаха синхронизации и смещает кривую BER вниз.
%     (б) С эквалайзером (сплошные линии): линейный фазовый набег внутри/около ЦП компенсируется по пилотам.
tags13_noEq = cell(1, numel(Mlist));
tags13_eq   = cell(1, numel(Mlist));
cfg13No = cfg; cfg13Eq = cfg;
if ~strcmpi(getenv('LAB3_QUICK'), '1')
    cfg13No.nFrames   = 300;
    cfg13No.minErrors = 120000;
    cfg13No.maxFrames = 3200;
    cfg13Eq.nFrames   = 100;
end

for iM = 1:numel(Mlist)
    % Без эквалайзера
    scenNo = scen_set(make_scen(), 'deltaMode', 'uniform', 'timeSync', 'cp', ...
                      'M', Mlist(iM), 'useEq', false);
    tagNo = sprintf('f13_noEq_M%d', Mlist(iM));
    outNo = ber_sweep(cfg13No, scenNo, cfg.snrDb, sprintf('без экв., M = %d', Mlist(iM)));
    S = add_curve(S, tagNo, outNo); tags13_noEq{iM} = tagNo;

    % С эквалайзером
    scenEq = scen_set(make_scen(), 'deltaMode', 'uniform', 'timeSync', 'cp', ...
                      'M', Mlist(iM), 'useEq', true);
    tagEq = sprintf('f13_eq_M%d', Mlist(iM));
    outEq = ber_sweep(cfg13Eq, scenEq, cfg.snrDb, sprintf('с экв., M = %d', Mlist(iM)));
    S = add_curve(S, tagEq, outEq); tags13_eq{iM} = tagEq;
end

h = new_figure(cfg, 'BER против M', 'Влияние когерентного накопления M и эквалайзера');
ax = axes(h); hold(ax, 'on'); grid(ax, 'on');
mkList = 'os^d';
for iM = 1:numel(Mlist)
    st = struct('col', cols(iM,:), 'mk', mkList(iM), 'ls', '-', 'open', false);
    plot_ber_curve(ax, get_curve(S, tags13_eq{iM}), st, sprintf('M = %d (с экв.)', Mlist(iM)));
end
plot_theory(ax, cfg);
for iM = 1:numel(Mlist)
    st = struct('col', cols(iM,:), 'mk', mkList(iM), 'ls', '--', 'open', true);
    plot_ber_curve(ax, get_curve(S, tags13_noEq{iM}), st, sprintf('M = %d (без экв.)', Mlist(iM)));
end
finish_ber_axes(ax, cfg, 'Влияние M и эквалайзера при синхронизации по ЦП');
legend(ax, 'Location', 'southwest', 'NumColumns', 2, 'FontSize', 8.5);
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
epsList = [0.01 0.1 0.3 0.5];
tags21 = cell(1, 2*numel(epsList));
for iE = 1:numel(epsList)
    for iEq = 1:2
        scen = scen_set(make_scen(), 'eps', epsList(iE), 'deltaMode', 'uniform', ...
                        'timeSync', 'genie', 'useEq', (iEq == 1), 'eqMode', 'cfo21');
        tag = sprintf('f21_e%d_eq%d', iE, iEq == 1);
        out = ber_sweep(cfg, scen, cfg.snrDb, sprintf('\\epsilon = %.2g, %s', ...
              epsList(iE), ternary(iEq == 1, 'с экв.', 'без экв.')));
        S = add_curve(S, tag, out); tags21{(iE-1)*2 + iEq} = tag;
    end
end

h = new_figure(cfg, 'BER при расстройке по частоте', 'Влияние дробной частотной расстройки');
ax = axes(h); hold(ax, 'on'); grid(ax, 'on');
colsE = [0 0.45 0.74; 0.30 0.75 0.93; 0.93 0.69 0.13; 0.64 0.08 0.18];
for iE = 1:numel(epsList)
    for iEq = 1:2
        st = struct('col', colsE(iE,:), 'mk', 'os^d', 'ls', ternary(iEq == 1, '-', '--'), ...
                    'open', iEq == 2);
        st.mk = st.mk(iE);
        plot_ber_curve(ax, get_curve(S, tags21{(iE-1)*2 + iEq}), st, ...
            sprintf('\\epsilon = %.2g, %s', epsList(iE), ...
            ternary(iEq == 1, 'с эквалайзером', 'без эквалайзера')));
    end
end
plot_theory(ax, cfg);
finish_ber_axes(ax, cfg, 'Без частотной синхронизации');
save_figure(cfg, h, 'fig_2_1_cfo');

% 2.2 Влияние M на частотную синхронизацию по ЦП (работа вблизи границы +-0.5)
cfg22 = cfg;
if ~strcmpi(getenv('LAB3_QUICK'), '1')
    cfg22.nFrames   = cfg.nFramesFine;
    cfg22.minErrors = 120000;
    cfg22.maxFrames = 2500;
end
eps22 = 0.496;
tags22 = cell(1, numel(Mlist));
for iM = 1:numel(Mlist)
    scen = scen_set(make_scen(), 'epsMode', 'fixed', 'eps', eps22, ...
                    'deltaMode', 'uniform', 'timeSync', 'genie', 'freqSync', 'cp', ...
                    'M', Mlist(iM), 'useEq', true);
    tag = sprintf('f22_M%d', Mlist(iM));
    out = ber_sweep(cfg22, scen, cfg.snrDb, sprintf('ЦП-ЧС, M = %d', Mlist(iM)));
    S = add_curve(S, tag, out); tags22{iM} = tag;
end

h = new_figure(cfg, 'BER частотной синхронизации по ЦП', 'Частотная синхронизация по ЦП');
ax = axes(h); hold(ax, 'on'); grid(ax, 'on');
mkList = 'os^d';
for iM = 1:numel(Mlist)
    st = struct('col', cols(iM,:), 'mk', mkList(iM), 'ls', '-', 'open', false);
    plot_ber_curve(ax, get_curve(S, tags22{iM}), st, sprintf('M = %d', Mlist(iM)));
end
cfgRef = cfg;  cfgRef.nFrames = min(cfg.nFrames, 20);
outNo22 = ber_sweep(cfgRef, scen_set(make_scen(), 'epsMode', 'fixed', ...
    'eps', eps22, 'deltaMode', 'uniform', 'timeSync', 'genie', ...
    'freqSync', 'none', 'useEq', true), cfg.snrDb, 'без синхр. частоты');
S = add_curve(S, 'f22_none', outNo22);
plot_ber_curve(ax, get_curve(S, 'f22_none'), ...
    struct('col', [0.5 0.5 0.5], 'mk', 'x', 'ls', ':', 'open', true), 'без синхронизации частоты');
plot_theory(ax, cfg);
finish_ber_axes(ax, cfg, sprintf('Частотная синхронизация по ЦП (\\epsilon = %.3f)', eps22));
save_figure(cfg, h, 'fig_2_2_cp_freq_M');

% 2.3 Влияние длины циклического префикса Ncp
%     При малых Ncp (4, 16) оценка частоты обладает большей дисперсией (СКО ~ 1/sqrt(Ncp)),
%     что у границы +-0.5 чаще приводит к срывам/перескокам и повышает BER.
NcpList = [4 16 64 256];
eps23 = 0.490;
tags23 = cell(1, numel(NcpList));

for iL = 1:numel(NcpList)
    cfgL = cfg;
    if ~strcmpi(getenv('LAB3_QUICK'), '1')
        cfgL.nFrames   = 240;
        cfgL.minErrors = 160000;
        cfgL.maxFrames = 2500;
    end
    cfgL.Ncp  = NcpList(iL);
    cfgL.Nsym = cfgL.Nfft + cfgL.Ncp;
    scen = scen_set(make_scen(), 'epsMode', 'fixed', 'eps', eps23, ...
                    'deltaMode', 'uniform', 'timeSync', 'genie', 'freqSync', 'cp', ...
                    'M', 1, 'useEq', true);
    tag = sprintf('f23_Ncp%d', NcpList(iL));
    out = ber_sweep(cfgL, scen, cfg.snrDb, sprintf('N_{CP} = %d', NcpList(iL)));
    S = add_curve(S, tag, out); tags23{iL} = tag;
end

h = new_figure(cfg, 'BER против длины ЦП', 'Влияние длины циклического префикса');
ax = axes(h); hold(ax, 'on'); grid(ax, 'on');
colsL = [0.85 0.33 0.10; 0.93 0.69 0.13; 0.47 0.67 0.19; 0 0.45 0.74];
mkList = 'os^d';
for iL = 1:numel(NcpList)
    st = struct('col', colsL(iL,:), 'mk', mkList(iL), 'ls', '-', 'open', false);
    plot_ber_curve(ax, get_curve(S, tags23{iL}), st, ...
        sprintf('N_{CP} = %d (%.1f%%)', NcpList(iL), 100*NcpList(iL)/cfg.Nfft));
end
plot_theory(ax, cfg);
finish_ber_axes(ax, cfg, sprintf('Частотная синхронизация по ЦП (\\epsilon = %.3f, M = 1)', eps23));
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
%LAB3_CFG  Единая точка настройки параметров моделирования.
    cfg.Nfft      = 2048;                 % размер БПФ (число поднесущих)
    cfg.Ncp       = 256;                  % длина ЦП (1/8 от Nfft)
    cfg.Nguard    = 256;                  % защитные полосы по краям (нулевые поднесущие)
    cfg.pilotStep = 4;                    % пилот на каждой 4-й поднесущей
    cfg.nDataSym  = 14;                   % OFDM-символов данных в кадре
    cfg.evalSym   = 3:12;                 % из них используются для оценки BER
    cfg.qamBits   = 4;                    % КАМ-16 -> 4 бита на символ
    cfg.snrDb     = 0:2:30;               % Eb/N0, дБ
    cfg.snrRef    = 14;                   % опорное Eb/N0 для сводной таблицы
    cfg.snrDemo   = 0;                    % Eb/N0 для демонстрации преамбулы и Рис. 0
    cfg.snrMetric12 = 10;                 % Eb/N0 для АКФ ЦП в п. 1.2 (как в эталоне)
    cfg.nFrames   = 40;                   % число кадров на точку BER (базовое)
    cfg.nFramesFine = 200;                % число кадров для тонких эффектов в п. 1.3, 2.2, 2.3
    cfg.nAvgMetric = 80;                  % усреднение метрики Минна по кадрам
    cfg.nAvgMetric12 = 50;                % усреднение АКФ по ЦП в п. 1.2
    cfg.snrDelta  = 12;                   % Eb/N0 для развёртки по delta
    cfg.nFramesDelta = 30;                % кадров на точку развёртки по delta
    cfg.seed      = 20241003;
    cfg.figureResolution = 170;

    % быстрый отладочный режим
    if strcmpi(getenv('LAB3_QUICK'), '1')
        cfg.nFrames = 4; cfg.nFramesFine = 4; cfg.snrDb = 0:6:12;
        cfg.nAvgMetric = 4; cfg.nAvgMetric12 = 4;
        cfg.nFramesDelta = 4;
    end

    cfg.Nsym  = cfg.Nfft + cfg.Ncp;
    cfg.padLen = cfg.Nfft + 1024 + 2*cfg.Nsym; % запас нулей от выхода за границы

    % --- Частотный план: центрированная сетка (-N/2 ... N/2 - 1) -----
    % DC-поднесущая (0 Гц) расположена в центре спектра (индекс Nfft/2 + 1)
    cfg.dcIdx = cfg.Nfft/2 + 1;

    % Защитные полосы (нули) слева (1..Nguard) и справа (Nfft-Nguard+1..Nfft)
    activeIdx = (cfg.Nguard + 1) : (cfg.Nfft - cfg.Nguard);

    % Пилоты (КАМ-4) с шагом pilotStep внутри активной полосы
    pilotCand = activeIdx(mod(activeIdx, cfg.pilotStep) == 0).';
    cfg.pilotIdx = setdiff(pilotCand, cfg.dcIdx);

    % Поднесущие данных (КАМ-16)
    cfg.dataIdx  = setdiff(activeIdx(:), union(cfg.pilotIdx, cfg.dcIdx));
    cfg.Nused    = numel(cfg.dataIdx) + numel(cfg.pilotIdx);

    % --- Известные последовательности (детерминированные) -------------
    % Пилоты: детерминированные известные символы КАМ-4 единичной мощности с плоским спектром
    pStream = RandStream('mt19937ar', 'Seed', 12345);
    pBits = randi(pStream, [0 1], 2, numel(cfg.pilotIdx));
    cfg.pilotSyms    = qammod(pBits, 4, 'InputType', 'bit', 'UnitAveragePower', true).';
    cfg.preambleTime = minn_preamble_time(cfg);                  % [A A -A -A]
    cfg.preambleX    = fftshift(fft(cfg.preambleTime));

    % --- Каталоги и видимость графиков --------------------------------
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
    X(:, 1) = cfg.preambleX;                       % преамбула (центрированный спектр)

    % Полезные биты: [4*numel(cfg.dataIdx), cfg.nDataSym]
    bits = randi([0 1], cfg.qamBits * numel(cfg.dataIdx), cfg.nDataSym);
    syms = qammod(bits, 16, 'InputType', 'bit', 'UnitAveragePower', true);
    X(cfg.dataIdx, 2:end)  = syms;
    X(cfg.pilotIdx, 2:end) = repmat(cfg.pilotSyms, 1, cfg.nDataSym);

    % Стандартное ОБПФ (без *sqrt(N)): средняя мощность отсчёта данных = Nused / Nfft^2
    x = ifft(ifftshift(X, 1), N, 1);
    xcp = [x(N-Ncp+1:N, :); x];                    % добавление циклического префикса
    frame = xcp(:);
    tx.bits = bits;
    tx.X = X;
end

% ------------------------------------------------------------------------
function [rx, ch] = lab_channel(cfg, frame, delta, epsNorm, ebnoDb)
%LAB_CHANNEL  Канал: временной сдвиг, частотный сдвиг, АБГШ.
%   Соглашение о знаке delta:
%     -Ncp <= delta < 0 - окно БПФ внутри ЦП: МСИ нет, только фазовый наклон;
%     delta < -Ncp      - окно выходит за пределы ЦП: появляется МСИ;
%     delta > 0         - окно запаздывает: МСИ из следующего символа.
    N = cfg.Nfft;  P = cfg.padLen;
    nTail = P + abs(delta) + 2*cfg.Nsym;
    rx = zeros(numel(frame) + P + nTail, 1);
    start = P + 1 - delta;                          % первый отсчёт кадра
    rx(start:start+numel(frame)-1) = frame;

    % частотный сдвиг: r[n] = s[n] * exp(j*2*pi*eps*n/N)
    n = (0:numel(rx)-1).';
    if epsNorm ~= 0
        rx = rx .* exp(1j*2*pi*epsNorm*n/N);
    end

    % АБГШ: формула с семинара с учетом единичной мощности созвездия (UnitAveragePower) и 1/Nfft
    Ps = mean(abs(qammod(0:15, 16, 'UnitAveragePower', true)).^2) / cfg.qamBits / N;
    N0 = Ps / (10^(ebnoDb / 10));
    sigma2 = N0;
    noise = sqrt(sigma2 / 2) * (randn(size(rx)) + 1j * randn(size(rx)));
    rx = rx + noise;

    ch.delta = delta; ch.eps = epsNorm; ch.sigma2 = sigma2;
end

% ------------------------------------------------------------------------
function [dHat, info] = lab_sync(cfg, scen, rx, delta)
%LAB_SYNC  Временная и частотная синхронизация приёмника.
%   dHat - оценка положения первого отсчёта ЦП первого символа кадра.
    dNom = cfg.padLen + 1;                          % ожидаемое начало кадра
    info.metric = []; info.axis = [];

    switch scen.timeSync
        case 'none'                                 % синхронизация не выполняется
            dHat = dNom;
        case 'genie'                                % точное положение (идеал)
            dHat = dNom - delta;
        case 'cp'                                   % слепая оценка по |АКФ| ЦП (как в timeSync.m)
            [dHat, info.metric, info.axis] = cp_metric(rx, cfg, scen.M, dNom, 'abs');
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

    bitsRef = tx.bits(:, cfg.evalSym);
    nBit = numel(bitsRef);

    if any(idx(:) < 1) || any(idx(:) > numel(rxc))
        % Катастрофический срыв синхронизации (окно за пределами массива)
        nErr = nBit;
        return;
    end

    % Приемное БПФ (без /sqrt(Nfft)) с переводом спектра в центрированный вид
    Y = fftshift(fft(rxc(idx), [], 1), 1);                       % N x nEval

    if scen.useEq
        if isfield(scen, 'eqMode') && strcmp(scen.eqMode, 'cfo21')
            H = lab_channel_estimate_21(cfg, Y);
        else
            H = lab_channel_estimate(cfg, Y);
        end
        Y = Y ./ H;
    end

    bitsHat = qamdemod(Y(cfg.dataIdx, :), 16, 'OutputType', 'bit', 'UnitAveragePower', true);
    nErr = sum(bitsHat(:) ~= bitsRef(:));
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
%LAB_CHANNEL_ESTIMATE  Оценка канала по пилотам (МНК, наклон ФЧХ и фаза символа).
%   1) МНК-оценка на пилотных поднесущих:
    Hp = Y(cfg.pilotIdx, :) ./ cfg.pilotSyms;        % [numPilots x nSym]
    pos = cfg.pilotIdx(:);
    grid = (1:cfg.Nfft).';

%   2) Оценка наклона ФЧХ (компенсация временного сдвига delta внутри ЦП):
%      Между соседними пилотами фаза набегает как exp(-j*2*pi*delta*k/N).
%      Оцениваем наклон slope, чтобы устранить сокращение хорды при интерполяции:
    posDiff = diff(pos);
    valid = (posDiff == cfg.pilotStep);             % исключаем разрыв вокруг DC
    dph = angle(Hp(2:end, :) .* conj(Hp(1:end-1, :)));
    slope = median(dph(valid, :) / cfg.pilotStep, 'all');

%   3) Устранение наклона фазы и усреднение формы канала по символам кадра:
    Hp_flat = Hp .* exp(-1j * slope * pos);
    Hm = mean(Hp_flat, 2);

%   4) Линейная интерполяция выпрямленной формы и возврат наклона ФЧХ:
    H_flat = interp1(pos, Hm, grid, 'linear', 'extrap');
    Hshape = H_flat .* exp(1j * slope * grid);

%   5) Оценка общей фазы каждого символа (CPE - Common Phase Error) для ЧС:
    cpe = mean(Hp ./ (Hshape(cfg.pilotIdx) * ones(1, size(Hp, 2))), 1);
    cpe = cpe ./ max(abs(cpe), realmin);             % сохраняем единичную амплитуду

%   6) Результирующая оценка канала для каждого символа:
    H = Hshape * cpe;
    mag = abs(H);
    bad = mag < 0.1;
    H(bad) = 0.1 * exp(1j * angle(H(bad)));
end

% ------------------------------------------------------------------------
function H = lab_channel_estimate_21(cfg, Y)
%LAB_CHANNEL_ESTIMATE_21  Локальная посимвольная оценка для сценария 2.1 (без ЧС).
%   Отслеживает фазовый набег каждого OFDM-символа при расстройке eps,
%   предотвращая взаимное погашение пилотов при усреднении по символам.
    Hp = Y(cfg.pilotIdx, :) ./ cfg.pilotSyms;        % [numPilots x nSym]
    pos = cfg.pilotIdx(:);
    grid = (1:cfg.Nfft).';

    posDiff = diff(pos);
    valid = (posDiff == cfg.pilotStep);             % исключаем разрыв вокруг DC
    dph = angle(Hp(2:end, :) .* conj(Hp(1:end-1, :)));
    slope = median(dph(valid, :) / cfg.pilotStep, 1); % [1 x nSym] посимвольно

    Hp_flat = Hp .* exp(-1j * (pos * slope));
    H_flat = interp1(pos, Hp_flat, grid, 'linear', 'extrap');
    H = H_flat .* exp(1j * (grid * slope));
    mag = abs(H);
    bad = mag < 0.1;
    H(bad) = 0.1 * exp(1j * angle(H(bad)));
end

% ------------------------------------------------------------------------
function [dHat, metric, axisD] = cp_metric(rx, cfg, M, dNom, metricMode)
%CP_METRIC  Слепая метрика временной синхронизации по ЦП с когерентным накоплением.
%   metricMode = 'abs' (по умолчанию) использует модуль |АКФ| (как в timeSync.m п. 3);
%   metricMode = 'rho' использует нормировку 2*|Aacc| / Pacc (для Рис. 0).
    if nargin < 5 || isempty(metricMode)
        metricMode = 'abs';
    end
    N = cfg.Nfft;  L = cfg.Ncp;  Nsym = cfg.Nsym;
    r = rx(:);
    r2 = abs(r).^2;

    c  = conj(r(1:end-N)) .* r(1+N:end);            % попарные произведения
    cs = [0; cumsum(c)];                            % префиксные суммы

    nc   = numel(c);
    dmax = nc - L + 1;
    dd   = (1:dmax).';
    A    = cs(dd+L) - cs(dd);                       % корреляция окна ЦП и хвоста

    if strcmpi(metricMode, 'rho')
        dValid = dmax - (M-1)*Nsym;
        if dValid < 1; error('cp_metric:short', 'Сигнал короче окна накопления'); end
        Aacc = A(1:dValid);
        for m = 1:M-1
            s = m*Nsym;
            Aacc = Aacc + A(1+s : dValid+s);
        end
        e  = [0; cumsum(r2)];
        Pn = (e(dd+L) - e(dd)) + (e(dd+N+L) - e(dd+N));
        Pacc = Pn(1:dValid);
        for m = 1:M-1
            s = m*Nsym;
            Pacc = Pacc + Pn(1+s : dValid+s);
        end
        P_thresh = 0.25 * max(Pacc);
        metric = 2 * abs(Aacc) ./ max(Pacc, P_thresh);
    else
        % Для 'abs' (п. 1.3) накапливаем по M случайным информационным OFDM-символам (m = 1..M),
        % чтобы при M = 1 не захватывалась детерминированная преамбула после нулевой паузы.
        dValid = dmax - M*Nsym;
        if dValid < 1; error('cp_metric:short', 'Сигнал короче окна накопления'); end
        Aacc = A(1+Nsym : dValid+Nsym);
        for m = 2:M
            s = m*Nsym;
            Aacc = Aacc + A(1+s : dValid+s);
        end
        metric = abs(Aacc) / M;
    end

    % Слепое обзорное окно по всей области возможного прихода кадра
    lo = max(1, dNom - floor(N/2) - 64);
    hi = min(dValid, dNom + ceil(N/2) + 64);
    [~, k] = max(metric(lo:hi));
    dHat   = lo + k - 1;                            % максимум корреляции
    axisD  = (lo:hi).' - dNom;
    metric = metric(lo:hi);
end

% ------------------------------------------------------------------------
function [dHat, metric, axisD] = minn_metric(rx, cfg, dNom)
%MINN_METRIC  Двухполовинная метрика Минна по преамбуле [A A -A -A] (из timeSync.m).
%   P = P1 + P2, Metric = |P|^2 / (E1 * E2) с порогом защиты от нулевой паузы.
    N = cfg.Nfft;  L = N/4;
    r = rx(:);

    c  = conj(r(1:end-L)) .* r(1+L:end);
    cs = [0; cumsum(c)];
    e  = [0; cumsum(abs(r).^2)];

    nc   = numel(c);
    dmax = nc - L + 1;
    dd   = (1:dmax).';
    P0   = cs(dd+L) - cs(dd);                       % корреляция соседних блоков длины L

    dValid = numel(r) - 4*L + 1;
    if dValid < 1; error('minn_metric:short', 'Сигнал короче преамбулы'); end

    P1 = P0(1:dValid);
    P2 = P0(1+2*L : dValid+2*L);
    P  = P1 + P2;

    E1 = e((1:dValid) + 2*L) - e(1:dValid);
    E2 = e((1:dValid) + 4*L) - e((1:dValid) + 2*L);

    lo = max(1, dNom + cfg.Ncp - floor(N/2) - 64);
    hi = min(dValid, dNom + cfg.Ncp + ceil(N/2) + 64);

    eWin = E2(lo:hi);
    eAct = median(eWin(eWin > 0.5 * max(eWin)));
    eFloor = 0.85 * eAct;
    metric = abs(P).^2 ./ (max(E1, eFloor) .* max(E2, eFloor));

    [~, k] = max(metric(lo:hi));
    dHat   = max(1, lo + k - 1 - cfg.Ncp);          % оценка начала кадра
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
%CFO_FROM_MINN  Оценка расстройки по преамбуле Минна (диапазон |eps|<2, как в freqSync.m).
    N = cfg.Nfft;  L = N/4;
    r = rx(:);
    i0 = dHat + cfg.Ncp;                            % начало тела преамбулы [A A -A -A]
    P1 = sum(conj(r(i0 + (0:L-1)))       .* r(i0 + L   + (0:L-1)));
    P2 = sum(conj(r(i0 + 2*L + (0:L-1))) .* r(i0 + 3*L + (0:L-1)));
    P = P1 + P2;
    epsHat = (2 / pi) * angle(P);
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
function y = minn_preamble_time(cfg)
%MINN_PREAMBLE_TIME  Преамбула [A A -A -A] с согласованной с данными мощностью.
    N = cfg.Nfft;  L = N/4;
    stream = RandStream('mt19937ar', 'Seed', 54321);
    b = randi(stream, [0 1], 2, L);
    q = qammod(b, 4, 'InputType', 'bit', 'UnitAveragePower', true).';
    A = ifft(q(:), L) * sqrt(L);                       % блок единичной мощности
    y = [A; A; -A; -A];
    y = y / sqrt(mean(abs(y).^2));                     % единичная средняя мощность
    y = y * sqrt(cfg.Nused) / N;                       % мощность равна мощности OFDM-символов
end

% ------------------------------------------------------------------------
function [acf, lags] = cyclic_acf(y)
%CYCLIC_ACF  Циклическая автокорреляционная функция (через спектр мощности).
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
%FINDPEAKS_SIMPLE  Локальные максимумы АКФ.
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
function [mAvg, axisD] = metric_shape(cfg, ~, M, ebnoDb, nAvg)
%METRIC_SHAPE  Усреднённая комплексная АКФ по ЦП (без деления на энергию и максимум).
%   Повторяет Пункт 2 эталона timeSync.m: ACF усредняется по M символам и по
%   nAvg реализациям до взятия модуля abs(ACF_avg).
    N = cfg.Nfft;  L = cfg.Ncp;  Nsym = cfg.Nsym;
    % Центрируем на информационном символе (символ 2: dNom = padLen + 1 + Nsym)
    dNom = cfg.padLen + 1 + cfg.Nsym;
    A_sum = []; axisD = [];
    for i = 1:nAvg
        txFrame = lab_frame(cfg);
        rx = lab_channel(cfg, txFrame, 0, 0, ebnoDb);
        r = rx(:);
        c  = conj(r(1:end-N)) .* r(1+N:end);
        cs = [0; cumsum(c)];
        nc = numel(c);
        dmax = nc - L + 1;
        dd = (1:dmax).';
        A = cs(dd+L) - cs(dd);
        dValid = dmax - (M-1)*Nsym;
        Aacc = A(1:dValid);
        for m = 1:M-1
            s = m*Nsym;
            Aacc = Aacc + A(1+s : dValid+s);
        end
        Aacc = Aacc / M;                            % усреднение по M символам
        lo = max(1, dNom - floor(1.4*L));
        hi = min(dValid, dNom + floor(1.4*L));
        if isempty(A_sum)
            A_sum = zeros(hi - lo + 1, 1);
            axisD = (lo:hi).' - dNom;
        end
        A_sum = A_sum + Aacc(lo:hi);
    end
    mAvg = abs(A_sum / nAvg);
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
function ber = ber_point(cfg, scen, ebnoDb, nFrames)
%BER_POINT  Расчёт BER при фиксированном Eb/N0.
    rng(cfg.seed);
    nErr = 0; nBit = 0;
    for nf = 1:nFrames
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
%BER_SWEEP  Расчёт кривой помехоустойчивости (BER от Eb/N0) с опциональным
%   адаптивным набором кадров по числу битовых ошибок.
    rng(cfg.seed);
    out.snr = snrDb(:).';
    out.ber = nan(size(snrDb));
    out.err = zeros(size(snrDb));
    out.bits = zeros(size(snrDb));
    out.label = label;

    hasMinErr = isfield(cfg, 'minErrors') && cfg.minErrors > 0;
    minF = cfg.nFrames;
    maxF = cfg.nFrames;
    minE = 0;
    if hasMinErr && isfield(cfg, 'maxFrames') && cfg.maxFrames > minF
        maxF = cfg.maxFrames;
        minE = cfg.minErrors;
    end

    for i = 1:numel(snrDb)
        nErr = 0; nBit = 0; nDone = 0;
        stopHighSnr = (i > 1) && ~isnan(out.ber(i-1)) && (out.ber(i-1) < 1.5e-5);
        frameCap = ternary(stopHighSnr, minF, maxF);
        while nDone < frameCap && (nDone < minF || (hasMinErr && nErr < minE))
            nDone = nDone + 1;
            [e, n] = frame_trial(cfg, scen, snrDb(i));
            nErr = nErr + e;  nBit = nBit + n;
        end
        out.err(i) = nErr;  out.bits(i) = nBit;
        if nBit > 0; out.ber(i) = nErr/nBit; end
        fprintf('  %-40s Eb/N0=%4.0f дБ: BER=%.3e (%4d/%6d, кадров %4d)\n', ...
                label, snrDb(i), out.ber(i), nErr, nBit, nDone);
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
%PLOT_THEORY  Теоретическая кривая помехоустойчивости КАМ-16 в АБГШ (berawgn).
    snr = linspace(cfg.snrDb(1), cfg.snrDb(end), 301);
    berTh = berawgn(snr, 'qam', 16);
    plot(ax, snr, berTh, 'k-.', 'LineWidth', 1.3, 'DisplayName', 'теория КАМ-16 (berawgn)');
end

% ------------------------------------------------------------------------
function finish_ber_axes(ax, cfg, ttl)
%FINISH_BER_AXES  Оформление осей BER (логарифмический масштаб по BER).
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
%SAVE_FIGURE  Сохранение рисунка в PNG без встроенной панели инструментов осей.
    axAll = findall(h, 'Type', 'axes');
    for i = 1:numel(axAll)
        try
            if isprop(axAll(i), 'Toolbar') && ~isempty(axAll(i).Toolbar)
                axAll(i).Toolbar.Visible = 'off';
            end
        catch
        end
    end
    drawnow;
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
    fid = fopen(fullfile(cfg.outDir, 'Lab3_selfcheck.txt'), 'w', 'n', 'UTF-8');
    frame = lab_frame(cfg);
    d0 = cfg.padLen + 1;
    N = cfg.Nfft; Ncp = cfg.Ncp; Nsym = cfg.Nsym;
    say = @(varargin) fprintf(fid, varargin{:});

    say('длина кадра: %d отсчётов (%.2f символа)\n', numel(frame), numel(frame)/Nsym);

    % 1) ЦП - точная копия хвоста символа (и проверка согласованности мощности)
    for s = 1:3
        off = (s-1)*Nsym;
        a = frame(off + (1:Ncp));
        b = frame(off + N + (1:Ncp));
        rho = abs(sum(conj(a).*b)) / sqrt(sum(abs(a).^2)*sum(abs(b).^2));
        say('символ %d: rho(ЦП, хвост) = %.6f, мощность = %.6e\n', ...
            s, rho, mean(abs(frame(off + (1:Nsym))).^2));
    end

    % 2) модулятор/демодулятор КАМ-16
    bt = randi([0 1], 4, 1000);
    syms = qammod(bt, 16, 'InputType', 'bit', 'UnitAveragePower', true);
    bh = qamdemod(syms, 16, 'OutputType', 'bit', 'UnitAveragePower', true);
    say('КАМ-16: BER(mod->demod) = %g (должно быть 0)\n', mean(bt(:) ~= bh(:)));

    % 3) метрика по ЦП при точном положении и без шума
    rx = lab_channel(cfg, frame, 0, 0, 100);
    [dHat, metric, ax] = cp_metric(rx, cfg, 1, d0, 'rho');
    say('метрика ЦП: dHat = %d (истина %d), rho(0) = %.4f\n', ...
        dHat, d0, metric(ax == 0));
    [dHatM, metM] = minn_metric(rx, cfg, d0);
    say('метрика Минна: dHat = %d (истина %d), max M = %.4f\n', ...
        dHatM, d0, max(metM));

    % 4) нормировка шума
    Ps = mean(abs(qammod(0:15, 16, 'UnitAveragePower', true)).^2) / cfg.qamBits / N;
    say('Eb/N0 = 10 дБ: Ps = %.6e, N0 = %.6e\n', Ps, Ps / 10);
    fclose(fid);
    fprintf('Самопроверка записана в %s\n', fullfile(cfg.outDir, 'Lab3_selfcheck.txt'));
end
