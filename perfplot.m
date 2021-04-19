data = readmatrix('data.csv');
N = data(:,1);
perf_naive = data(:, 2);
perf_strassen = data(:, 3);
perf_hybrid = data(:, 4);

plot1 = loglog(N, perf_naive, 'LineWidth', 2);
plot1.Color(4) = 0.7;
hold on
plot2 = plot(N, perf_strassen, 'LineWidth', 2);
plot2.Color(4) = 0.7;
plot3 = plot(N, perf_hybrid, 'LineWidth', 2);
plot3.Color(4) = 0.7;
title('Benchmarking')
xlabel('Input size (NxN matrices)')
ylabel('Runtime (ms)')
legend('Naive', 'Strassen', 'Hybrid')

% x=0:0.1:10;
% y1=x.^3;
% y2=x.^2.8704;
% y3=x.^2.3729;
% plot(x,y1)
% title('Asymptotic complexity comparison')
% xlabel('Input size (NxN matrices)')
% ylabel('Runtime')
% hold on
% plot(x,y2)
% plot(x,y3)
% legend('O(x^3) - Naive', 'O(x^2.8704) - Strassen', 'O(x^2.3729) - Coppersmith-Winograd')