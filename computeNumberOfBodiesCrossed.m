function nbScatterers = computeNumberOfBodiesCrossed(TXpos, RXpos, people2DPositionA, people2DPositionB)

% Incident ray
%people2DPositionA = [people2DPositionA(1,1,1,1), people2DPositionA(1,1,1,2)];
nbScatterers1 = computeBC(TXpos, people2DPositionA, people2DPositionB);

%Reflected ray
nbScatterers2 = computeBC(people2DPositionA, RXpos, people2DPositionB);

nbScatterers = nbScatterers1 + nbScatterers2;

end

function nBC = computeBC(startPoint, endPoint, peoplePoints)
nBC = 0;
distance_finder = @(point1, point2) sqrt(sum((point1 - point2).^2));
lineDistance = distance_finder(startPoint, endPoint);
for i = 1:size(peoplePoints, 1)
    currentPerson = peoplePoints(i, :);
    start2point = distance_finder(startPoint, currentPerson);
    point2end = distance_finder(currentPerson, endPoint);
    radius = 0.15;
    AC = start2point;
    BC = point2end;
    AB=lineDistance;
    u = endPoint-currentPerson;
    v = startPoint - currentPerson;
    CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    teta = real(acosd(CosTheta));
    alpha = asin(BC*sind(teta)/AB);
    d = abs(AC*sin(alpha));
    if d<=radius
      nBC = nBC + (radius-d)/radius;
    end
    
end
end