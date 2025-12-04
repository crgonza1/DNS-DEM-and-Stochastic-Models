function [x_new pdf_new] = sum_pdfs(x1,pdf1,x2,pdf2)

dx      = x1(2) - x1(1);
x_new   = linspace(min(x1)+min(x2),max(x1)+max(x2),length(x1)+length(x2)-1);
pdf_new = conv(pdf1,pdf2);
pdf_new = pdf_new/trapz(x_new,pdf_new);

return