% Create a movie from the phi data
phis_im_100_002 = importdata('phis_im_100_0.02.dat');
phis_re_100_002 = importdata('phis_re_100_0.02.dat');

phis_im_100_005 = importdata('phis_im_100_0.05.dat');
phis_re_100_005 = importdata('phis_re_100_0.05.dat');

phis_im_200_002 = importdata('phis_im_200_0.02.dat');
phis_re_200_002 = importdata('phis_re_200_0.02.dat');

phis_im_200_005 = importdata('phis_im_200_0.05.dat');
phis_re_200_005 = importdata('phis_re_200_0.05.dat');

L = 2;
xmin = -1.0;
N = size(phis_re_200_005,1) - 1;
x = (L/N) .* (0:N) + xmin;

% Make a plot for each state vector in phi
for k = 1:size(phis_re_200_005,2)
    plot(x, phis_re_200_005(:,k) );
    xlabel('x', 'FontSize', 16);
    ylabel('Re(\phi)', 'FontSize', 16);
    title('Wave Packet propagating in Free Space with k = 200, \sigma = 0.05', 'FontSize', 18);
    axis([-1, 1, -1, 1]);
    M_re_200_005(k) = getframe;
end
movie2gif(M_re_200_005, 'M_re_200_005.gif')

for k = 1:size(phis_re_200_002,2)
    plot(x, phis_re_200_002(:,k) );
    xlabel('x', 'FontSize', 16);
    ylabel('Re(\phi)', 'FontSize', 16);
    title('Wave Packet propagating in Free Space with k = 200, \sigma = 0.02', 'FontSize', 18);
    axis([-1, 1, -1, 1]);
    M_re_200_002(k) = getframe;
end
movie2gif(M_re_200_002, 'M_re_200_002.gif')
movie2avi(M_re_200_002, 'M_re_200_002.avi')


for k = 1:size(phis_re_100_005,2)
    plot(x, phis_re_100_005(:,k) );
    xlabel('x', 'FontSize', 16);
    ylabel('Re(\phi)', 'FontSize', 16);
    title('Wave Packet propagating in Free Space with k = 100, \sigma = 0.05', 'FontSize', 18);
    axis([-1, 1, -1, 1]);
    M_re_100_005(k) = getframe;
end
movie2gif(M_re_100_005, 'M_re_100_005.gif')
movie2avi(M_re_100_005, 'M_re_100_005.avi')


for k = 1:size(phis_re_100_002,2)
    plot(x, phis_re_100_002(:,k) );
    xlabel('x', 'FontSize', 16);
    ylabel('Re(\phi)', 'FontSize', 16);
    title('Wave Packet propagating in Free Space with k = 100, \sigma = 0.02', 'FontSize', 18);
    axis([-1, 1, -1, 1]);
    M_re_100_002(k) = getframe;
end
movie2gif(M_re_100_002, 'M_re_100_002.gif')
movie2avi(M_re_100_002, 'M_re_100_002.avi')

phis_rev_im_200_005 = importdata('phis_rev_im_200_0.05.dat');
phis_rev_re_200_005 = importdata('phis_rev_re_200_0.05.dat');
% Make a plot for each state vector in phi
for k = 1:size(phis_rev_re_200_005,2)
plot(x, phis_rev_re_200_005(:,k) );
xlabel('x', 'FontSize', 16);
ylabel('Re(\phi)', 'FontSize', 16);
title('Wave Packet propagating in Free Space before Reversal with k = 200, \sigma = 0.05', 'FontSize', 18);
axis([-1, 1, -1, 1]);
M_rev_re_200_005(k) = getframe;
end
movie2gif(M_rev_re_200_005, 'M_rev_re_200_005.gif')
movie2avi(M_rev_re_200_005, 'M_rev_re_200_005.avi')
