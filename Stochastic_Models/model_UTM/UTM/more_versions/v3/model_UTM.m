function [t,xp_mean,xp_var,xp_skew,xp_kur,btc1,btc2,btc3,btc4,btc5] = model_UTM(np,max_sim_time,dt,xi,pdf_xi,diameter)


% Compute BTCs:
[btc1,btc2,btc3,btc4,btc5] = get_btc(max_sim_time,dt,xi,pdf_xi,diameter);
%toc
%plot(btc1(:,1),btc1(:,2),btc2(:,1),btc2(:,2),btc3(:,1),btc3(:,2),btc4(:,1),btc4(:,2),btc5(:,1),btc5(:,2))
%keyboard


% Compute the statistics of the particles motion
[t,xp_mean,xp_var,xp_skew,xp_kur] = get_stats(np,max_sim_time,dt,xi,pdf_xi);


return
