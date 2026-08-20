%% ================================================================
%  Moment-rate distributions (Total vs Locked) 
%  Inputs expected in workspace:
%    rates        : [n_patches x 1]  long-term slip rate      [mm/yr]
%    creep_rates  : [n_realizations x n_patches] creep rate   [mm/yr]
%    locked_index : [n_realizations x n_patches] locking (0-1)
%    patch_stuff.area_faces : [n_patches x 1] patch area      [km^2]
%
%  Physics:  Mdot0 = mu * A * (v_plate - v_creep)
%            Mw    = (2/3)*(log10(M0[N.m]) - 9.1)   (Hanks & Kanamori 1979)
% =================================================================

% ---------- constants ----------
mu      = 3.0e10;                 % shear modulus [Pa]
mm2m    = 1e-3;                   % mm -> m
yr2sec  = 365.25*24*3600;         % yr -> s
area_m2 = patch_stuff.area_faces(:).' * 1e6;   % km^2 -> m^2, row [1 x n_patches]

% ---------- back-slip (slip deficit) per realization ----------
bslip = rates(:).' - creep_rates;              % [n_real x n_patch], mm/yr
bslip_m_per_s = bslip * mm2m / yr2sec;         % mm/yr -> m/s

% ---------- moment-rate deficit per patch, per realization ----------
moment_deficit_all = mu .* area_m2 .* bslip_m_per_s;   % [N.m/s]
moment_deficit_all = moment_deficit_all * yr2sec / 1e19;   % -> 10^19 N.m/yr

% ---------- sum over patches -> one value per realization ----------
total_moment_rate  = sum(moment_deficit_all, 2);
locked_moment_rate = sum(moment_deficit_all .* locked_index, 2);

% =================================================================
%  1) Histograms
% =================================================================
figure;
histogram(total_moment_rate, 'BinMethod','auto', ...
          'FaceColor',[0 0.7 0.6], 'EdgeColor','k');
xlabel('Total Moment Rate (10^{19} N·m/yr)');
ylabel('Frequency');
title('Total Moment Rate Histogram');
grid on; box on;

figure;
histogram(locked_moment_rate, 'BinMethod','auto', ...
          'FaceColor',[0.6 0.4 1], 'EdgeColor','k');
xlabel('Locked Moment Rate (10^{19} N·m/yr)');
ylabel('Frequency');
title('Locked Patch Moment Rate Histogram');
grid on; box on;

% =================================================================
%  2) KDE plot with Mw rows (100 / 300 / 1000 yr) + colored heading
% =================================================================
fig = figure('Units','centimeters','Position',[2 2 26 15]);
ax  = axes(fig, 'Position', [0.08 0.11 0.90 0.60]);   % leave headroom
hold(ax,'on');

x = linspace(2, 22, 400);

% --- Total ---
f_tot  = ksdensity(total_moment_rate, x);
hTotal = fill(ax, [x fliplr(x)], [f_tot zeros(size(f_tot))], [0 0.70 0.60], ...
     'FaceAlpha',0.4, 'EdgeColor','none');
plot(ax, x, f_tot, 'Color',[0 0.50 0.40], 'LineWidth',2, ...
     'HandleVisibility','off');

% --- Locked ---
f_lock  = ksdensity(locked_moment_rate, x);
hLocked = fill(ax, [x fliplr(x)], [f_lock zeros(size(f_lock))], [0.60 0.40 1], ...
     'FaceAlpha',0.4, 'EdgeColor','none');
plot(ax, x, f_lock, 'Color',[0.40 0 1], 'LineWidth',2, ...
     'HandleVisibility','off');

xlim(ax, [x(1) x(end)]);
ylim(ax, [0 max([f_tot f_lock])*1.05]);

xlabel(ax, 'Moment rate (10^{19} N m yr^{-1})', 'FontSize',18);
ylabel(ax, 'Probability density',               'FontSize',18);
legend([hTotal hLocked], {'Total moment rate','Locked moment rate'}, ...
       'Location','northwest', 'FontSize',16, 'AutoUpdate','off');
set(ax, 'FontSize',14, 'Box','on');

% -----------------------------------------------------------------
%  Mw tick rows — drawn as text/line primitives (robust, no
%  undocumented properties, works in all recent releases)
% -----------------------------------------------------------------
years   = [100, 300, 1000];
colors  = {[0.85 0.47 0.00], [0.12 0.56 1.00], [0.70 0.13 0.13]};
rowGap  = [0.06, 0.16, 0.26];        % row offsets, fraction of axis height

drawnow;                              % finalize ticks/limits
xt = ax.XTick;  xt = xt(xt > 0);      % Mw undefined for <= 0
yl = ylim(ax);  H = diff(yl);  yTop = yl(2);
xL = xlim(ax);

% spine + tick marks for the FIRST row only (matches the Python plot)
ySpine  = yTop + rowGap(1)*H;
tickLen = 0.015*H;
line(ax, xL, [ySpine ySpine], 'Color','k', 'LineWidth',0.8, ...
     'Clipping','off', 'HandleVisibility','off');
for xv = xt
    line(ax, [xv xv], [ySpine ySpine+tickLen], 'Color','k', ...
         'Clipping','off', 'HandleVisibility','off');
end

% Mw labels, one row per recurrence interval
for k = 1:numel(years)
    Mw   = (2/3) * (log10(xt * 1e19 * years(k)) - 9.1);
    yRow = yTop + rowGap(k)*H + 1.3*tickLen;
    for j = 1:numel(xt)
        text(ax, xt(j), yRow, sprintf('%.2f', Mw(j)), ...
             'Color', colors{k}, 'FontSize',13, ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', 'Clipping','off');
    end
end

% -----------------------------------------------------------------
%  Centred heading with colored year numbers (single TeX string)
% -----------------------------------------------------------------
yHead = yTop + (rowGap(end) + 0.13)*H;
headTxt = ['Magnitude ranges based on moment rates (' ...
    '{\color[rgb]{0.85,0.47,0.00}100}, ' ...
    '{\color[rgb]{0.12,0.56,1.00}300}, ' ...
    '{\color[rgb]{0.70,0.13,0.13}1000} yr)'];
text(ax, mean(xL), yHead, headTxt, ...
     'Interpreter','tex', 'FontSize',18, ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', 'Clipping','off');