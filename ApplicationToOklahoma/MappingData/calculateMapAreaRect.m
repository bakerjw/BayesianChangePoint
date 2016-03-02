function area = calculateMapAreaRect(localVect)
% This function calculates the area within a local region
% defined by localVect - [siteLat, siteLong, deltaLat, deltaLong].
% deltaLat and deltaLong define the +- ranges around siteLAt and siteLong.


siteLat = localVect(1);
siteLong = localVect(2);
dLat = localVect(3);
dLong = localVect(4);
    
lat1 = siteLat - dLat; lat2 = siteLat + dLat;
long1 = siteLong - dLong; long2 = siteLong + dLong;

earthEllipsoid = referenceSphere('earth','km');
area = areaquad(lat1, long1, lat2, long2, earthEllipsoid);