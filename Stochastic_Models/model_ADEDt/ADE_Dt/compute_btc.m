function [hist_BTCs] = compute_btc(C,x,diameter,hist_BTCs,flag_xtrm)

last_x = x(end-flag_xtrm);

% Compute pdfs of the time crossing dxups(1,2,...5). 
if last_x>=0.1*diameter
    d     = 0.1*diameter;
    [a,b] = search_indexs(x,d);
    c     = (d-x(a))/(x(b)-x(a));
    hist_time_dxups1 = C(:,a) + (C(:,b)-C(:,a))*c;
else
    hist_time_dxups1 = C(:,1)*NaN;
end

if last_x>=1*diameter 
    d     = 1*diameter;
    [a,b] = search_indexs(x,d);
    c     = (d-x(a))/(x(b)-x(a));
    hist_time_dxups2 = C(:,a) + (C(:,b)-C(:,a))*c;
else
    hist_time_dxups2 = C(:,1)*NaN;
end

if last_x>=10*diameter
    d     = 10*diameter;
    [a,b] = search_indexs(x,d);
    c     = (d-x(a))/(x(b)-x(a));
    hist_time_dxups3 = C(:,a) + (C(:,b)-C(:,a))*c;    
else
    hist_time_dxups3 = C(:,1)*NaN;
end

if last_x>=100*diameter
    d     = 100*diameter;
    [a,b] = search_indexs(x,d);
    c     = (d-x(a))/(x(b)-x(a));
    hist_time_dxups4 = C(:,a) + (C(:,b)-C(:,a))*c;    
else
    hist_time_dxups4 = C(:,1)*NaN;
end

if last_x>=1000*diameter
    d     = 1000*diameter;
    [a,b] = search_indexs(x,d);
    c     = (d-x(a))/(x(b)-x(a));
    hist_time_dxups5 = C(:,a) + (C(:,b)-C(:,a))*c;    
else
    hist_time_dxups5 = C(:,1)*NaN;
end

hist_BTCs = [hist_BTCs; [hist_time_dxups1 hist_time_dxups2 hist_time_dxups3 hist_time_dxups4 hist_time_dxups5]];

return