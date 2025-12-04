function [btc1,btc2,btc3,btc4,btc5] = get_btc(max_sim_time,dt,xi,pdf_xi,diameter)

% Do you want to fix the amount of elements of the pdf? (useful when it has a lot of elements)
flag_fix = 1;                                                              

if flag_fix == 1
    nu_elm = 100;
    if nu_elm < length(xi)        
        xi_fix = linspace(min(xi),max(xi),nu_elm);
        pdf_xi = interp1(xi,pdf_xi,xi_fix);
        
        xi     = xi_fix(:);
        pdf_xi = pdf_xi(:)/trapz(xi,pdf_xi);
    end
end


% Get an uniform spatial distribution
dx = min(diff(xi));
if dx~=max(diff(xi))
    xi_aux = min(xi):dx:max(xi);
    pdf_xi = interp1(xi,pdf_xi,xi_aux);
    xi     = xi_aux;    
end
pdf_xi = pdf_xi/trapz(xi,pdf_xi);


% Compute pdf_xi for each t
t          = [dt:dt:max_sim_time]';                                        % Set the time domain
xi_old     = 0;
pdf_old    = 1;
zero_thres = 10^-15;
BTC        = zeros(length(t),5);
timer1     = cputime;
threshold = 10^-100;

for i=1:length(t)                                                          % Loop over t
    
    % Remove small values of the current pdf
    if min(pdf_old)<threshold
        [xi_old pdf_old] = remove_near_zero(xi_old,pdf_old,threshold);     % Remove small values of the pdf
    end
    
    % Compute the sum of pdf_0 + pdf_t
    [xi_new pdf_new] = sum_pdfs(xi,pdf_xi,xi_old,pdf_old);
    
    
    for j=1:5
        dxups = (10^(j-2))*diameter;
        if min(xi_new)<=dxups && dxups <=max(xi_new)
            BTC(i,j) = interp1(xi_new,pdf_new, dxups);
        end        
    end
                
    xi_old   = xi_new;
    pdf_old  = pdf_new;
    
    % Display message
    if mod(i,floor(length(t)/100))==0
        timer2 = cputime;
        disp([num2str(round(100*i/length(t))),'% completed of the BTCs (current proccess took ',num2str(timer2-timer1),' s).'])
        timer1 = timer2;
    end
end


% Normalizate and arrange the BTCs
for i=1:5
    if trapz(t,BTC(:,i))==0
        a(i) = 1;
    else
        a(i) = trapz(t,BTC(:,i));
    end
end

btc1(:,1) = t(:); btc1(:,2) = BTC(:,1)/a(1);
btc2(:,1) = t(:); btc2(:,2) = BTC(:,2)/a(2);
btc3(:,1) = t(:); btc3(:,2) = BTC(:,3)/a(3);
btc4(:,1) = t(:); btc4(:,2) = BTC(:,4)/a(4);
btc5(:,1) = t(:); btc5(:,2) = BTC(:,5)/a(5);

return