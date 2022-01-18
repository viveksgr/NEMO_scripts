function shaded_plot(x,y_mean, y_std,k)
% Plots y_mean(1xN) against x(1xN) with y_std (1xN or 1x1) for shading around the mean
% line y_mean.
if nargin<4
    k = 'r'; % Default color is red
end
curve1 = y_mean + y_std;
curve2 = y_mean - y_std;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, k, 'FaceAlpha',0.2,'EdgeColor','None');
hold on;
plot(x, y_mean, 'Color', k, 'LineWidth', 1);
end