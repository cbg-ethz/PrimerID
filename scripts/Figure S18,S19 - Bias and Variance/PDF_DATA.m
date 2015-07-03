function logPr = PDF_DATA(Data_raw, Data_pID, Data_HC, theta)
logDirPDF = @(p, mu, s) gammaln(s) - sum(gammaln(mu*s)) + sum((mu*s-1).*log(p));

mu_raw = theta(1:5);
mu_pID = theta(6:10);
mu_HC  = theta(11:15);

 s_raw = theta(16);
 s_pID = theta(17);
 s_HC  = theta(18);

IND = [1, 2, 3, 4, 5, 6];

logPr = 0;
for i = IND
    % Raw frequencies
    logPr = logPr + logDirPDF(Data_raw(i,:)', mu_raw, s_raw);
    
    % pID frequencies
    logPr = logPr + logDirPDF(Data_pID(i,:)', mu_pID, s_pID);
    
    %  HC frequencies
    logPr = logPr + logDirPDF(Data_HC(i,:)' , mu_HC , s_HC );
end
end

