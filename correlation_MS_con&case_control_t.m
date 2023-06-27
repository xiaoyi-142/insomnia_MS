a = 0;
b = 0;
c = 0;
d = 0;

xvalues = meanMS_con;
yvalues = meantstat;

nregs = length(xvalues);

for ind = 1:nregs
    xval = xvalues(ind);
    yval = yvalues(ind);
    if (xval < 0) && (yval > 0) % a is top left quadrant
        a = a + 1;
    elseif (xval > 0) && (yval > 0) % b is top right quadrant
        b = b + 1;
    elseif (xval < 0) && (yval < 0) % c is bottom left quadrant
        c = c + 1;
    elseif (xval > 0) && (yval < 0) % d is bottom right quadrant
        d = d + 1;
    end
end

a_percentage = a / nregs; % percentage of scatter points in the top left quadrant
b_percentage = b / nregs; % percentage of scatter points in the top right quadrant
c_percentage = c / nregs; % percentage of scatter points in the bottom left quadrant
d_percentage = d / nregs; % percentage of scatter points in the bottom right quadrant

a_percentage = round(a_percentage * 100, 2); % percentage of scatter points in the top left quadrant
b_percentage = round(b_percentage * 100, 2); % percentage of scatter points in the top right quadrant
c_percentage = round(c_percentage * 100, 2); % percentage of scatter points in the bottom left quadrant
d_percentage = round(d_percentage * 100, 2); % percentage of scatter points in the bottom right quadrant

% 绘制散点图和回归线
figure
scatter(meanMS_con, meantstat, 'x', 'k')
refline(0, 0)
hold on
plot([0 0], ylim)
plot([0 0], ylim, '-b')

% 绘制回归线
coefficients = polyfit(meanMS_con, meantstat, 1);
regression_line = polyval(coefficients, meanMS_con);
plot(meanMS_con, regression_line, 'b')

% 添加a、b、c、d的标记
text(0.65 * min(meanMS_con), 0.55 * max(meantstat), [sprintf('%.2f', a_percentage) '%'], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontName','Arial','FontSize',20)
text(0.65 * max(meanMS_con), 0.55 * max(meantstat), [sprintf('%.2f', b_percentage) '%'], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontName','Arial','FontSize',20)
text(0.65 * min(meanMS_con), 0.55 * min(meantstat), [sprintf('%.2f', c_percentage) '%'], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontName','Arial','FontSize',20)
text(0.65 * max(meanMS_con), 0.55 * min(meantstat), [sprintf('%.2f', d_percentage) '%'], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontName','Arial','FontSize',20)

correlation = corrcoef(meanMS_con, meantstat);
correlation = correlation(1, 2);
correlation = round(correlation, 2); % 保留两位小数

% 添加相关系数文本到参考线旁边
text(0.65 * max(meanMS_con), 0.15 * max(meantstat), ['r =  ' sprintf('%.2f', correlation)], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontName','Arial','FontSize',20)

xlabel('Mean control MS', 'FontName','Arial','FontSize',25)
ylabel('Case-control t-statistic', 'FontName','Arial','FontSize',25)

% 调整坐标轴的字号大小
set(gca, 'FontSize', 20)

hold off

% 设置保存图像的 DPI
dpi = 300; % 设置为所需的 DPI

% 保存图片到指定位置，并指定 DPI
print('results/correlation_image.png', '-dpng', ['-r' num2str(dpi)]);
