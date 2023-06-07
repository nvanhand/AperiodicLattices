%user inputs
Xrange = [-20,20];
Yrange = [-20,20];

sigma_input = .2;

Nx = 5;
Ny = 5;

Ntotal = Nx*Ny;
%calculate total length of the bounding box
Xlength = Xrange(2)-Xrange(1);
Ylength = Yrange(2)-Yrange(1);

%calculate cell size
Cx = Xlength/Nx;
Cy = Ylength/Ny;

sigma = sigma_input*(Cx*0.5+Cy*0.5);
%create a loop to store the variables
pointsmatrix = zeros(Ntotal,2); %blank matrix to hold points
pointsmatrix_rnd = zeros(Ntotal,2);

%create a uniform periodic point set
counter = 1;
for i=1:Nx
    for j=1:Ny
        xi = i*Cx;
        yj = j*Cy;
        pointsmatrix(counter,1) = xi;
        pointsmatrix(counter,2) = yj;
        counter = counter + 1;
    end
end

%add perturbation to the existing point set
for i=1:Ntotal
    point_i = pointsmatrix(i,:);
    
    %create the random vector
    xr = normrnd(0,sigma);
    yr = normrnd(0,sigma);
    
    perturbationvector = [xr,yr];
    
    %apply the perturbation to the point
    point_r = point_i + perturbationvector;
    
    pointsmatrix_rnd(i,:) = point_r;
    
end
pointsmatrix_rnd
%plot the points
figure
hold on
scatter(pointsmatrix(:,1),pointsmatrix(:,2));
scatter(pointsmatrix_rnd(:,1),pointsmatrix_rnd(:,2),'*');
