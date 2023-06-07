%the below code just creates a square polygon with two holes. Any polygon
%input you want to try should be called "newpoly"
testmode = true;
if testmode==true
Center = [0,0];
Width = 5;
Height = 10;
P1 = Center + [-Width/2,Height/2];
P2 = Center + [Width/2,Height/2];
P3 = Center + [-Width/2,-Height/2];
P4 = Center + [Width/2,-Height/2];
Boundarymatrix = [P1;P2;P4;P3];
polygon_in = polyshape(Boundarymatrix(:,1),Boundarymatrix(:,2)); 
 smallbox = [[0,0];[1,0];[1,1];[0,1]];
 smallbox2 = smallbox + [1.2,1.2];
 smallpoly = polyshape(smallbox(:,1),smallbox(:,2));
 smallpoly2 = polyshape(smallbox2(:,1),smallbox2(:,2));
 newpoly = subtract(subtract(polygon_in,smallpoly),smallpoly2);
 thickness = 5; %defines thickness of the shape
end
filename = 'test_1.stl';
 
stlholder = poly2_stl(newpoly,thickness);
stlwrite(stlholder,filename)

function [testtri] = poly2_stl(polyin,thickness)
%this function should be able to take in a 2D polyshape with holes, and
%export a .STL file which can then be 3D printed.


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%user input


%%get the holes of the polyshape and store into a cell array
polyholes = holes(polyin); %holes of the polygon
polybound = rmholes(polyin); %outline of the polygon
allpolys = [polyholes;polybound];

for i=1:length(allpolys)
    %get .STL boundary of allpolys(i)
    inlooptri = boundary2stlwall2(allpolys(i),thickness);
    %add .STL boundary to the blob
    if i==1
        triBlob = inlooptri;
    else
        triBlob = combineTris2(triBlob,inlooptri);
    end
end
%%create a 2D triangulation of newpoly
tri = triangulation(polyin);



points = tri.Points;
Triangles = tri.ConnectivityList;

points2_z = zeros(length(points),1); %contains z values for the first face/shell of the polygon .stl
points3_z = zeros(length(points),1)+thickness; %contains z valuse for the second face/shell of the polygon .stl

%initialize blank matrices to hold the point surfaces of the shape
points3 = [points,points3_z];
points2 = [points,points2_z];

triface1 = triangulation(Triangles,points3);

%flip the triangles to make the bottom face normals not inside out
Trianglesreversed = Triangles;
Trianglesreversed(:,2) = Triangles(:,3);
Trianglesreversed(:,3) = Triangles(:,2);

triface2 = triangulation(Trianglesreversed,points2);
trifaces = combineTris2(triface1,triface2);

testtri = combineTris2(triBlob,trifaces);


end

function [triout] = boundary2stlwall2(polyin,t)
%this function takes in a single closed polyshape and turns it into a
%triangulation and point set which can be used to create a .STL file from
%the boundary. t is the desired thickness of the wall.

%%%%%%%%%%%%%%%%%testing stuff
%pointvals = [[1.1,2.5];[2.6,3.2];[3.3,4.7];[2.9,-5.9];[-5.9,-5.9]]; %points to create the polyshape
%randomshape = polyshape(pointvals(:,1),pointvals(:,2)); %creating a dummy polyshape
%t = 1.2;

[X,Y] = boundary(polyin); %takes the polyshape and separates the edge points into x and y matrices 
pointvals = [X,Y];

%each polyshape will have a triangulation matrix listing where the
%triangles are connected. This triangulation will contain two triangles for
%every link in the polyshape. For every link there is one point
tripoints = zeros(length(pointvals)*2,3);
triangles = zeros(length(pointvals)*2,3);
%triangles = [0,0,0];
counter = 1;
for i=1:length(pointvals)
    if i<length(pointvals)
        %[pointvals(i,:) ; pointvals(i+1,:)]
        pA = pointvals(i,:);
        pB = pointvals(i+1,:);
        
    else
        %[pointvals(i,:) ; pointvals(1,:)]
        pA = pointvals(i,:);
        pB = pointvals(1,:);
    end
        [points,tris] = link2tri2(pA,pB,t);
        
        %add these new points to the tripoints matrix
        tripoints(counter,:) = points(1,:);
        tripoints(counter+1,:) = points(2,:);
        tripoints(counter+2,:) = points(3,:);
        tripoints(counter+3,:) = points(4,:);
        counter = counter + 4;
        
        %add these tris to the triangles matrix
        
        tris = tris + (i-1)*4;
        triangles(i*2-1,:) = tris(1,:);
        triangles(i*2,:) = tris(2,:);
        
   
end

%create a triangulation of the now combined triangles and points
triout = triangulation(triangles,tripoints);

%write to a .stl file
%stlwrite(testtri,'testyboi.stl')
end

function [triout] = combineTris2(tri1,tri2)
%takes in two triangulations and combines them into one
points1 = tri1.Points;
points2 = tri2.Points;
Triangles1 = tri1.ConnectivityList;
Triangles2 = tri2.ConnectivityList;

pointscombined = [points1;points2];
trianglescombined = [Triangles1;Triangles2+length(points1)];

triout = triangulation(trianglescombined,pointscombined);
end

function [points,tris] = link2tri2(pA,pB,t)
%this function takes in two points and the desired thickness and generates
%two triangles and a set of points for it.
%pA = [1,2];
%pB = [2,3];
%t = 5; 

p1 = [pA 0];
p2 = [pB 0];
p3 = [pA t];
p4 = [pB t];


tris = [[1 3 2];[3 4 2]]; %change to 1 2 3 if shit gets fucked
points = [p1;p2;p3;p4];

end