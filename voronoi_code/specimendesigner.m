type_list = ['s','h','t'];
wall_list = [0.5,0.6,0.7];
sigma_list = [0,0.05,0.1,0.15,.2];
perturbation_type_list = ['seed','node'];
sigma_node_list = [0,0.05,0.1,0.2];


total = length(sigma_list)*length(perturbation_type_list),length(type_list),length(wall_list),length(type_list);

%points = createSeeds2([0,0],4,2,4,2,'t');
%scatter(points(:,1),points(:,2));




%user inputs (primary)
partnumber = 66;
strutthickness = .7;
sigma_nucleus = 0;
sigma_node = 0;
type = 'h';
voronoimaker5(strutthickness,sigma_nucleus,sigma_node,type,partnumber)

Samples_table = readtable('Sample_List.csv');
Samples_variables = Samples_table.Variables;
count = 0;

densitymatrix = zeros(length(Samples_variables),2)

for i=1:length(Samples_variables)

    partnumber = Samples_variables(i,1);
    strutthickness = Samples_variables(i,2);
    sigma_node = Samples_variables(i,3);
    sigma_nucleus = Samples_variables(i,4);

    density = voronoimaker5(strutthickness,sigma_nucleus,sigma_node,type,partnumber);
    count = count + 1;
    partnumber
    densitymatrix(i,:) = [partnumber,density];
end
csvwrite("specimen_rel_density.csv",densitymatrix);