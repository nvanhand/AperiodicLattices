strutthickness = 0.7;  %mm
sigma_nucleus = 0; 
sigma_node = 0.2; 
type = 't'; 

ending = ['V', string(sigma_nucleus), '_N', string(sigma_node), '_', type, '-'];

trials = 100;

area_data = zeros(trials+2,2);
area_data(:, 1) = 0:trials+1;

% Unperturbed 
baseshape = voronoimaker5(strutthickness, 0, 0,type, '');
%plotshape(finalshape, 'triangle_honeycomb.png');
area_data(2,2) = baseshape.area;

% Iterate versions of Perturbed
for n=2:(trials+2)
filename = join([ending string(n) '.png'], ''); 
[finalshape] = voronoimaker5(strutthickness,sigma_nucleus,sigma_node,type, n);
%plotshape(finalshape, filename)
area = finalshape.area;
area_data(n+1,2) = area;
end

writematrix(area_data, 'area.csv')

function plotshape(data, filename)
f = figure('visible','off');
plot(data)
saveas(f, filename)
close(f);
end
