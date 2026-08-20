

%faulttrace=fliplr([-120.62 36.06;-120.49 35.945;-120.37 35.84;-120.26 35.7]);
%origin = fliplr([-120.4337 35.9395 ]);
%faulttrace = llh2localxy(faulttrace', origin);

%faulttrace=[-16.8327   13.3702; -5.0942    0.6103; 5.7714  -11.0399; 18.2765 -24.6333];
%faulttrace=[-30.3702   28.0921; 18.2765  -24.6333];
%faulttrace=[-30.3352   28.0429; 14.52  -20.8];
faulttrace=[ -23.5686   20.6799; 21.2866  -28.1630];
%faulttrace=[-30.3352   28.0429; 21.2866  -28.1630];

%enter fault geometry
	%faults=[length, width, *depth, dip, strike(degrees), *north offset, *east offset]
	% *depth to top edge, north and east offsets refer to location of center of top edge
    length1=sqrt((faulttrace(1,1)-faulttrace(2,1))^2+(faulttrace(1,2)-faulttrace(2,2))^2);
	strike1=90-atan((faulttrace(2,2)-faulttrace(1,2))/(faulttrace(2,1)-faulttrace(1,1)))*180/pi;
    centerx1=(faulttrace(1,1)+faulttrace(2,1))/2;
    centery1=(faulttrace(1,2)+faulttrace(2,2))/2;
    
    %length1 = length1*2;
    
    dip=-90;
    W=15;
    faults(1,:)=[length1 W 0 dip strike1  centerx1 centery1];
    
    
pm=[];
pmss=[];
pmds=[];

for k=1:size(faults,1)
    %specify components of slip to be calculate ([strike-slip,dip-slip,opening]) -- e.g. [0 1 0] means dip slip only
dis_geom1  = [faults(k,:), [1 1 0]];
dis_geom = movefault(dis_geom1);  % move the fault so that the coordinates of the midpoint refer to the
											 % fault bottom as in Okada

pm1=patchfault(dis_geom(1,1:7),nhe,nve);

%% Create slip patches
for j=1:nhe
pm = [pm; pm1(1+nve*(j-1):nve*j,:)];
end %j

end
