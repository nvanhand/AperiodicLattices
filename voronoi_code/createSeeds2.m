function [Pointvals] = createSeeds2(Center,width,height,Nx,Ny,type,strutthickness)

%temporary parameters for testing%%%%%%%%%%%%%%

%example input
%createSeeds2([0,0],4,2,4,2,'s')


%width = 4;
%height = 2;
%Nx = 4;
%Ny = 2;
%type = 's'createSeeds2([0,0],4,2,4,2,'s')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P1 = Center + [width*0.5,height*0.5]; %Northeast
P2 = Center + [width*0.5,-height*0.5]; %Northwest
P3 = Center - [width*0.5,height*0.5]; %Southwest
P4 = Center + [-width*0.5,height*0.5]; %Southeast
Box = [P1;P2;P3;P4];
Boxpoly = polyshape(Box(:,1),Box(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this function takes in parameters, such as box width, box height, center
%location, number of points in the X direction, number of points in the Y
%direction, and the type of seed dispersal such as square grid, triangular
%grid, and hexagonal grid.

%different types of patterns
% 's' = "square"
% 't' = "triangle"
% 'h' = "hexagonal"

%the buffer tells you how many rows of points will be extended outside the
%box boundary.
bufferx = 3;
buffery = 3;
%factor to scale nodes which fall outside the box boundary. set to 1 if you don't want any nodepushing. 
%the effect is that it helps eliminate weird borders at the edge of the
%specimen. May have unexpected effects if your center isn't at [0,0];
nodepush = 1; 
test = 0;
counter = 0;

Nx_ext = Nx + bufferx*2;
Ny_ext = Ny + buffery*2;
Ntotal = Nx_ext*Ny_ext;


cellwidth = width/Nx;
cellheight = height/Ny;

Pointvals = zeros(Ntotal,2);
if type == 's'
    %generates square seed patter

    for i=0:Nx_ext-1
        for j=0:Ny_ext-1
            counter = counter + 1;
           
            Pointvals(counter,1) = i*cellwidth;
            Pointvals(counter,2) = j*cellheight;
        end
    end
    R = [-(width*0.5+cellwidth*bufferx-cellwidth*(0.5)),-(height*0.5+cellheight*buffery-cellheight*0.5)];
    for i=1:length(Pointvals)
        Pointvals(i,:) = Pointvals(i,:) + R;
    end
elseif type =='t'
    %generates hexagonal seed pattern

    cellheight = cellwidth*sqrt(3)*0.5; %correct cell height for hexagonal honeycomb
    for i=0:Nx_ext-1
        for j=0:Ny_ext-1
            counter = counter + 1;
           if mod(j,2)==0
               xshift = cellwidth*0.5;
           else
               xshift = 0;
           end
            Pointvals(counter,1) = i*cellwidth + xshift;
            Pointvals(counter,2) = j*cellheight;
            counter = counter + 1;
        end
        
    end
    R = [-(width*0.5 + cellwidth*bufferx - cellwidth*0.5),...
        -(height*0.5 + cellheight*buffery - cellheight + strutthickness*0.5)];
    for i=1:length(Pointvals)
        Pointvals(i,:) = Pointvals(i,:) + R;
    end
    Pointvals = unique(Pointvals,'rows');
    [Vertices,~] = voronoin(Pointvals);
    for i=1:length(Vertices)-1
        if abs(Vertices(i,1))>width*3
            Vertices(i,:) = [];
        elseif abs(Vertices(i,2))>width*3
            Vertices(i,:) = [];
        end
    end
    Pointvals = Vertices;

%apply nodepushing to points which fall outside of the box boundary
for i=1:length(Pointvals)
    if ~isinterior(Boxpoly,Pointvals(i,:))
        Pointvals(i,:) = Pointvals(i,:)*nodepush;
     
    end
end


elseif type == 'h'
    %generates hexagonal seed pattern
    Rcenter = [0,0];
    cellheight = cellwidth*sqrt(3)*0.5; %correct cell height for hexagonal honeycomb
    for i=0:Nx_ext-1
        for j=0:Ny_ext-1
            counter = counter + 1;
           if mod(j,2)==0
               xshift = cellwidth*0.5;
           else
               xshift = 0;
           end
            Pointvals(counter,1) = i*cellwidth + xshift;
            Pointvals(counter,2) = j*cellheight;
            counter = counter + 1;
        end
        
    end
    R = [-(width*0.5+cellwidth*bufferx-cellwidth*(0.5)),-(3.4631+height*0.5-cellheight*buffery-cellheight*0.5+cellheight*Ny*0.5)];
    for i=1:length(Pointvals)
        Pointvals(i,:) = Pointvals(i,:) + R;
    end
  
end