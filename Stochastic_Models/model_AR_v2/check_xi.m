% This function moves xi (if mean xi <=0) in order to match mean_xp

function [xi] = check_xi(xi,pdf_xi,mean_xp)

ve = trapz(xi,xi.*pdf_xi)/trapz(xi,pdf_xi);

if ve <= 0
	xi = xi - ve + mean_xp;
    disp('  ')
    disp('  ')
    disp('  ')
    disp('WARNING: THE MEAN TRAVEL STEP PER DT_UPS')
    disp('COMPUTED FROM XI AND PDF_XI IS NEGATIVE!.')
    disp(['NOW IT WAS CHANGED FROM ',num2str(ve),' TO ',num2str(mean_xp),'.'])
    disp('  ')
    disp('  ')
    disp('  ')
end

return