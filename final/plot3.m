% Create a movie from the phi data

phis_step_re_200_20000 = importdata('phis_step_re_200_20000.0.dat');
phis_step_re_200_10000 = importdata('phis_step_re_200_10000.0.dat');
phis_step_re_200_30000 = importdata('phis_step_re_200_30000.0.dat');

L = 2;
xmin = -1.0;
N = size(phis_step_im_100_10,1) - 1;
x = (L/N) .* (0:N) + xmin;

% Make a plot for each state vector in phi
for k = 1:size(phis_step_re_200_20000,2)
    plot(x, phis_step_re_200_20000(:,k) );
    xlabel('x', 'FontSize', 16);
    ylabel('Re(\phi)', 'FontSize', 16);
    title('Wave Packet propagating in Free Space with k = 200, \sigma = 0.05', 'FontSize', 18);
    axis([-1, 1, -1, 1]);
    M_step_re_200_20000(k) = getframe;
end
movie2gif(M_step_re_200_20000, 'M_step_re_200_20000.gif')
%movie2avi(M_step_re_200_20000, 'M_step_re_200_20000.avi')

for k = 1:size(phis_step_re_200_10000,2)
    plot(x, phis_step_re_200_10000(:,k) );
    xlabel('x', 'FontSize', 16);
    ylabel('Re(\phi)', 'FontSize', 16);
    title('Wave Packet propagating in Free Space with k = 200, \sigma = 0.05', 'FontSize', 18);
    axis([-1, 1, -1, 1]);
    M_step_re_200_10000(k) = getframe;
end
movie2gif(M_step_re_200_10000, 'M_step_re_200_10000.gif')
%movie2avi(M_step_re_200_10000, 'M_step_re_200_10000.avi')

for k = 1:size(phis_step_re_200_30000,2)
    plot(x, phis_step_re_200_30000(:,k) );
    xlabel('x', 'FontSize', 16);
    ylabel('Re(\phi)', 'FontSize', 16);
    title('Wave Packet propagating in Free Space with k = 200, \sigma = 0.05', 'FontSize', 18);
    axis([-1, 1, -1, 1]);
    M_step_re_200_30000(k) = getframe;
end
movie2gif(M_step_re_200_30000, 'M_step_re_200_30000.gif')
%movie2avi(M_step_re_200_30000, 'M_step_re_200_30000.avi')