function [finalshape] = voronoimaker5(strutthickness,sigma_nucleus,sigma_node,type,partnumber)
%Create filename based on input parameters
%strutthickness = 2;
%sigma_nucleus = 0.1;
%sigma_node = 0;
%type = 'h';
%example input: voronoimaker5(0.5,0.05,0,'h',"12")
%ending = ["wall",strutthickness,"mm","_sv",sigma_nucleus,"_sn",sigma_node,"_type-",type,".csv"];
ending = ['V', string(sigma_nucleus), '_N', string(sigma_node), '_', type, '-', partnumber, '.csv'];
%Create filename based on input parameters for the overall structure
filename = join(ending, '');

%user inputs(secondary)
filletradius = 0;
filletradius_eff = filletradius*2;
Center = [0,0];
Width = 80;
Height = 80;
Nx = 10;
Ny = 10;

generatenew = true;
writedata = true;

%change parameters for triangle mode
if type=='t'
    Nx = 10;
    Ny = 10;
    Height = Width*sqrt(3)/2-strutthickness*0.5;
elseif type=='h'
    Height = Height - (Width/Nx)*sqrt(3)/2;
    %delete this line later
    Height = 69.2820;
end

%correct for cell size
sigma_nucleus = sigma_nucleus*Width/Nx;
sigma_node = sigma_node*Width/Nx;

%define outer boundary shape 
outerbound = points2rec([-(Width-strutthickness)*0.5,0],[(Width-strutthickness)*0.5,0],Height-strutthickness);
negativespace = outerbound;


%create voronoi seeds
if generatenew
    Points = createSeeds2(Center,Width,Height,Nx,Ny,type,strutthickness);
end
Points = unique(Points,'rows');

%apply nucleus perturbation
if generatenew
    randomvectors = zeros(length(Points),2);
end
for i=1:length(Points)
    if generatenew
        randomvectors(i,:) = [normrnd(0,sigma_nucleus),normrnd(0,sigma_nucleus)];
    end
    Points(i,1) = Points(i,1) + randomvectors(i,1);
    Points(i,2) = Points(i,2) + randomvectors(i,2);
end
%create voronoi tesselation
[V,C] = voronoin(Points);

%apply node perturbation to the corners
randomvectors_node = zeros(length(V),2);
for i=1:length(V)
    if generatenew
    randomvectors_node(i,1) = normrnd(0,sigma_node);
    randomvectors_node(i,2) = normrnd(0,sigma_node);
    end
    
    V(i,1) = V(i,1) + randomvectors_node(i,1);
    V(i,2) = V(i,2) + randomvectors_node(i,2);
end

%create a polygon which contains circles at all the corners
for i=1:length(V)

    if (i==1)
        circlepolygons = circlePoly([V(i,1),V(i,2)],strutthickness*0.5);
    else
        inloopcircle =  circlePoly([V(i,1),V(i,2)],strutthickness*0.5);
        circlepolygons = union(circlepolygons,inloopcircle);
    end
    
    
end
count = 1;

%now convert the polygons into individual struts, and combine them into one
%single polyshape
strutconnections = zeros(length(V)*6,2); %empty matrix to hold strut connection values

for i=1:length(C)
    inloopoly = cell2mat(C(i)); %convert the polygon in index i of cell array C into a normal matrix
    %checkervalue = sum([isnan(inloopoly),isinf(inloopoly)]);
    %if checkervalue==0
        %inloopolyshape = polyshape(inloopoly(:,1),inloopoly(:,2));
        %polyarea = area(inloopolyshape);
    %end
    
  %if polyarea>threshold
      %check to make sure the polygon is big enough to plot
    for j=1:length(inloopoly)
        %loop through the poly and create a strut connection for each side
        if j==length(inloopoly)
           pindex1 = inloopoly(j);
           pindex2 = inloopoly(1);
        
        else
           pindex1 = inloopoly(j);
           pindex2 = inloopoly(j+1);
           
        end
           strutconnections(count,1) = pindex1;
           strutconnections(count,2) = pindex2;
           point1 = V(pindex1,:);
           point2 = V(pindex2,:);
           
            %make sure it only creates a strut if the strut falls within a
            %reasonable range (Width+Height)*3) in order to create the
            %strut
           if abs(V(pindex1,1))>Width*3
               
           elseif abs(V(pindex1,2))>(Width+Height)*3
               
           elseif abs(V(pindex2,1))>(Width+Height)*3
               
           elseif abs(V(pindex2,2))>(Width+Height)*2
               
           else
               %if it's not too big, create the recangle
               rectangle = points2rec(point1,point2,strutthickness);
               
               %now the individual rectangle strut can be subtracted from the
           %negative space polyshape
           negativespace = subtract(negativespace,rectangle);
           end       
            
    end
    
  %end %end of polygon area check
end
negativespace = subtract(negativespace,circlepolygons);
if ~filletradius==0
        negativespace = polybuffer(negativespace,-filletradius_eff);
        negativespace = polybuffer(negativespace,filletradius_eff);
end

outer = polybuffer(outerbound,strutthickness,'JointType','miter');

%subtract negative space to create the final shape
finalshape = subtract(outer,negativespace);

%plot(finalshape)
%disp(finalshape)
%csvwrite(filename, finalshape);
%density = area(finalshape)/area(outer);

%axis square
%hold on
%scatter(Points(:,1),Points(:,2));
%plot(negativespace);
%write the final shape to a .STL file
stlholder = poly2stl(finalshape,extrusiondepth);
%%%%only use this if you want to scale the final .STL file
stlholder = scaleTris(stlholder,(1/1000)); %scale back to mm
stlwrite(stlholder,filename)

end %end of whole fucking thing