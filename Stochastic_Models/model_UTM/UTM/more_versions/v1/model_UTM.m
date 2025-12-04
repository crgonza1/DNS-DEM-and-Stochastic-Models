function [t,x,xp_mean,xp_var,btc1,btc2,btc3,btc4,btc5,C] = model_UTM(np,max_sim_time,dt,xi,pdf_xi,diameter)

% Define the time discretization
t = 0:dt:max_sim_time;

% Generate random length steps for each dt
rand_xi      = randpdf(xi,pdf_xi,np,length(t)-1);

% Compute the particles position for each particle
aux          = cumsum(rand_xi,2);
positions    = [zeros(np,1) aux];

% Compute the particles mean position as a function of the time
xp_mean = mean(positions);

% Compute the variance of the particles displacement
for i=1:length(t)
    xp_var(i) = var(positions(:,i)-xp_mean(i));
end

% Generate the spatial discretization
Ex = trapz(xi,xi.*pdf_xi)/trapz(xi,pdf_xi); % expected travel step during dt
dx = Ex/20;
x  = min(positions(:,2)):dx:min(positions(:,end));

% Add the BTC points to the spatial discretization:
dxups = diameter*[0.1 1 10 100 1000];
for i=1:5
   k = find(x==dxups(i));      
   if isempty(k)
       flag1 = 0;
   else
       flag1 = 1;
   end
   
   if x(end)>=dxups(i)
      flag2 = 1;
   else
      flag2 = 0;
   end
   
   if flag1==0 && flag2==1
      x = sort([x dxups(i)]); 
   end
end

% Compute the travel time for each particle
for i=1:np    
   A = positions(i,:);     
   diffA = [diff(A) 0];
        
   for j=1:length(x)
       a = find((x(j)-A).*(A+diffA-x(j))>=0); a=min(a);
       if isempty(a)
           a=NaN;
           b=NaN;
       elseif a==length(A);
           a=length(A)-1;
           b=length(A);
       else
           b=a+1;
       end
              
       if isnan(a)
           travel_times(i,j) = NaN;           
       else
           travel_times(i,j) = dt*(a-1) + dt*(x(j)-A(a))/(A(b)-A(a));
       end
       
   end
end

% Definition of BTCs:
for i=1:5
    if dxups(i)<=x(end)
        k = find(x==dxups(i));
        Naux = histc(travel_times(:,k),t);
        if i==1
            btc1(:,1) = t;
            btc1(:,2) = Naux/trapz(t,Naux);            
        elseif i==2
            btc2(:,1) = t;
            btc2(:,2) = Naux/trapz(t,Naux);
        elseif i==3
            btc3(:,1) = t;
            btc3(:,2) = Naux/trapz(t,Naux);
        elseif i==4
            btc4(:,1) = t;
            btc4(:,2) = Naux/trapz(t,Naux);
        elseif i==5
            btc5(:,1) = t;
            btc5(:,2) = Naux/trapz(t,Naux);
        end                                    
    end   
end

% Compute particles concentration
for j=1:length(x)
    C(:,j) = histc(travel_times(:,j),t);
end

return