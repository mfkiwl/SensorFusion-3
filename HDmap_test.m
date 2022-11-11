% %% Test using Open Street Map
% name = 'openstreetmap';
% url = 'https://a.tile.openstreetmap.org/${z}/${x}/${y}.png';
% copyright = char(uint8(169));
% attribution = copyright + "OpenStreetMap contributors";
% addCustomBasemap(name,url,'Attribution',attribution)
% 
% %% 
% zoomLevel = 12;
% player = geoplayer(sol.full.opt.P_lla(1,1),sol.full.opt.P_lla(1,2),zoomLevel);
% player.Basemap = 'openstreetmap';
% plotRoute(player,sol.full.opt.P_lla(:,1),sol.full.opt.P_lla(:,2));

%% Load HD Map Data (NGII)
%
% HD Map Data obtained from: http://map.ngii.go.kr/ms/pblictn/preciseRoadMap.do
%
% cd D:/SJ_Dataset/HDMap/Map1/HDMap_UTMK_타원체고/
cd D:/SJ_Dataset/HDMap/Map2/SEC01_BRT_내부간선도로/HDMap_UTMK_타원체고/
% cd D:/SJ_Dataset/HDMap/Map2/SEC02_세종정부청사_주변/HDMap_UTMK_타원체고/

T1 = readgeotable("A1_NODE.shp");
T2 = readgeotable("A2_LINK.shp");

%% Process data
Node = table2struct(T1);
Link = table2struct(T2);
NodeID = zeros(length(Node),1);

for i=1:length(Node)
    Node(i).ID = str2num(Node(i).ID{1}(end-4:end));
    NodeID(i) = Node(i).ID;
end

for i=1:length(Link)
    Link(i).ID = str2num(Link(i).ID{1}(end-4:end));
    
    if isempty(Link(i).R_LinkID{1})
        Link(i).R_LinkID = -1;
    else
        Link(i).R_LinkID = str2num(Link(i).R_LinkID{1}(end-4:end));
    end
    
    if isempty(Link(i).L_LinkID{1})
        Link(i).L_LinkID = -1;
    else
        Link(i).L_LinkID = str2num(Link(i).L_LinkID{1}(end-4:end));
    end

    Link(i).FromNodeID = str2num(Link(i).FromNodeID{1}(end-4:end));
    Link(i).ToNodeID = str2num(Link(i).ToNodeID{1}(end-4:end));
end

%% Convert Points to [Lat,Lon] Format
% UTM-K to WGS84
x_b = 1e6; y_b = 2e6;
NodeLLA = zeros(length(Node),2);
for i=1:length(Node)
    Node_lla = enu2lla([Node(i).Shape.X - x_b,Node(i).Shape.Y - y_b,0],[38,127.5,0],'ellipsoid');
    Node(i).lat = Node_lla(1); Node(i).lon = Node_lla(2);
    NodeLLA(i,1) = Node_lla(1); NodeLLA(i,2) = Node_lla(2);
end

%% Plot Data using Geoplot!

figure(1); 
geoplot(NodeLLA(:,1),NodeLLA(:,2),'r+'); hold on; grid on;
for i=1:length(Link)
    idxF = find(NodeID == Link(i).FromNodeID);
    idxT = find(NodeID == Link(i).ToNodeID);
    NodeIdxs = [idxF, idxT];
    geoplot(NodeLLA(NodeIdxs,1),NodeLLA(NodeIdxs,2),'k-')
    % Cannot represent left/right turn intervals and other road types
end

figure(2);
mapshow(T1);
mapshow(T2); grid on; axis equal; 

% Conclusion: Converting UTM-K to WGS84 should not be performed.
% Other than straight lane links, curved lane links cannot be properly
% transferred.
% When comparing, convert vehicle trajectory and lane related coordinates
% to UTM-K for accurate HD-map lane plot