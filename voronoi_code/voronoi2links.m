function [links] = voronoi2links(P,C)
%this function takes in two arrays containing information to build a set of
%voronoi, P and C and turns them into a set on "links" as a 2x2xN array of
%point values. one pair of points is one link and that contains the
%information needed to connect a single strut when creating a voronoi based
%structure. Links with infinite values are deleted and points that are
%closer to eachother than the threshold are merged. P is the points matrix,
%and C is the cell array containing the polygons
len = length(C);
count = 0;
%first we need to figure out the max number of links needed
for i=1:len
    thispoly = cell2mat(C(i));
    numlines = length(thispoly)+1;
    count = numlines + count;
end
links = zeros(2,2,count);
%now we can fill the links matrix with actual point values
count = 1;
for i=1:len
    thispoly = cell2mat(C(i));
    for j=1:length(thispoly)
        if j<length(thispoly)
        pindex1 = thispoly(j);
        pindex2 = thispoly(j+1);
        else 
        pindex1 = thispoly(length(thispoly));
        pindex2 = thispoly(1);
        end
        p1 = P(pindex1,:);
        p2 = P(pindex2,:);
        count = count + 1;
        pair = [p1;p2];
        links(:,:,count)=pair;
    end
    
    
end
end