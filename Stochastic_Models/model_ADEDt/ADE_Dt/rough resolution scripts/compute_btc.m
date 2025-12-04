function [hist_BTCs] = compute_btc(C,x,diameter,hist_BTCs,flag_xtrm)

last_x = x(end-flag_xtrm);

% Compute pdfs of the time crossing dxups(1,2,...5). 
if last_x>=0.1*diameter
    hist_time_dxups1 = C(:,x==0.1*diameter);
else
    hist_time_dxups1 = C(:,1)*NaN;
end

if last_x>=1*diameter    
    hist_time_dxups2 = C(:,x==1*diameter);
else
    hist_time_dxups2 = C(:,1)*NaN;
end

if last_x>=10*diameter    
    hist_time_dxups3 = C(:,x==10*diameter);
else
    hist_time_dxups3 = C(:,1)*NaN;
end

if last_x>=100*diameter    
    hist_time_dxups4 = C(:,x==100*diameter);
else
    hist_time_dxups4 = C(:,1)*NaN;
end

if last_x>=1000*diameter    
    hist_time_dxups5 = C(:,x==1000*diameter);
else
    hist_time_dxups5 = C(:,1)*NaN;
end

hist_BTCs = [hist_BTCs; [hist_time_dxups1 hist_time_dxups2 hist_time_dxups3 hist_time_dxups4 hist_time_dxups5]];

return