pore_radius = 1.0;
pore_length = 5.0;
bulk_radius = 10.0;
bulk_length = 10.0;

dxReservoir = 0.1;
dxReservoirWall = 0.1;
dxChannelCenter = 0.1;
dxChannelWall = 0.1;

Point(1) = {0, 0, 0, dxReservoir};
Point(2) = {bulk_length, 0, 0, dxChannelCenter};
Point(3) = {bulk_length+pore_length, 0, 0, dxChannelCenter};
Point(4) = {bulk_length+pore_length, pore_radius, 0, dxChannelWall};
Point(5) = {bulk_length, pore_radius, 0, dxChannelWall};
Point(6) = {bulk_length, bulk_radius, 0, dxReservoirWall};
Point(7) = {0, bulk_radius, 0, dxReservoirWall};
Point(8) = {-bulk_length, bulk_radius, 0, dxReservoirWall};
Point(9) = {-bulk_length, pore_radius, 0, dxChannelWall};
Point(10) = {-bulk_length-pore_length, pore_radius, 0, dxChannelWall};
Point(11) = {-bulk_length-pore_length, 0, 0, dxChannelCenter};
Point(12) = {-bulk_length, 0, 0, dxChannelCenter};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};

Line Loop(13) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Plane Surface(13) = {13};
