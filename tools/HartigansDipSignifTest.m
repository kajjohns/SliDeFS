
function		[dip, p_value, xlow,xup]=HartigansDipSignifTest(xpdf,nboot,varargin)

%  function		[dip, p_value, xlow,xup]=HartigansDipSignifTest(xpdf,nboot)
%
% calculates Hartigan's DIP statistic and its significance for the empirical p.d.f  XPDF (vector of sample values)
% This routine calls the matlab routine 'HartigansDipTest' that actually calculates the DIP
% NBOOT is the user-supplied sample size of boot-strap
% Code by F. Mechler (27 August 2002)

% calculate the DIP statistic from the empirical pdf

plotid = 6;
plothist = 0;
j = 1;
while j < nargin -1
    if strncmpi(varargin{j},'plot',4)
        plothist = 1;
        if isnumber(varargin{j+1})
            plotid = varagin{j+1};
            j = j+1;
        else
            plotid = 1;
        end
    end
    j = j+1;
end
[dip,xlow,xup, ifault, gcm, lcm, mn, mj]=HartigansDipTest(xpdf);
N=length(xpdf);

% calculate a bootstrap sample of size NBOOT of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
boot_dip=[];
for i=1:nboot
   unifpdfboot=sort(unifrnd(0,1,1,N));
   [unif_dip]=HartigansDipTest(unifpdfboot);
   boot_dip=[boot_dip; unif_dip];
end;
boot_dip=sort(boot_dip);
p_value=sum(dip<boot_dip)/nboot;

% Plot Boot-strap sample and the DIP statistic of the empirical pdf
if plothist
    figure(plotid); clf;
    [hy,hx]=hist(boot_dip);
    bar(hx,hy,'k'); hold on;
    plot([dip dip],[0 max(hy)*1.1],'r:');
end