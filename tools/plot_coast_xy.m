function handles = plot_coast_xy(bbox,origin,color)
%bbox [minlon minlat;maxlon maxlat]


if nargin<3
    color = 'k';
end

S = shaperead('ne_10m_coastline','BoundingBox',bbox);

for k=1:length(S)
    
    X= S(k).X;
    Y=S(k).Y;
    xy = llh2local([X;Y],origin)';
    handles(k) = plot(xy(:,1),xy(:,2),'color',color);
    
end
