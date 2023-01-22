function W = threshold_bybin(W, lowbin)
%THRESHOLD_PROPORTIONAL     Proportional thresholding
%
%   W_thr = threshold_proportional(W, p);
%
%   This function "thresholds" the connectivity matrix by preserving a
%   proportion p (0<p<1) of the strongest weights. All other weights, and
%   all weights on the main diagonal (self-self connections) are set to 0.
%
%   Inputs: W,      weighted or binary connectivity matrix
%           p,      proportion of weights to preserve
%                       range:  p=1 (all weights preserved) to
%                               p=0 (no weights preserved)
%
%   Output: W_thr,  thresholded connectivity matrix
%
%
%   Mika Rubinov, U Cambridge,
%   Roan LaPlante, Martinos Center, MGH
%   Zitong Zhang, Penn Engineering

%   Modification history:
%   2010: Original (MR)
%   2012: Bug fix for symmetric matrices (RLP)
%   2015: Improved symmetricity test (ZZ)
%   2019: Mitsouko van Assche: adapted for windowed procedure 
%   This function "thresholds" the connectivity matrix by selecting
%   connections according to a threshold range . All other weights, and
%   all weights on the main diagonal (self-self connections) are set to 0.

n=size(W,1);                                %number of nodes
W(1:n+1:end)=0;                             %clear diagonal

if max(max(abs(W-W.'))) < 1e-10             %if symmetric matrix
    W=triu(W);                              %ensure symmetry is preserved
    ud=2;                                   %halve number of removed links
else
    ud=1;
end

ind=find(W);                                %find all links
nbl = round((numel(ind)/10));             %number of links in each bin
E=sortrows([ind W(ind)], 2);                %sort by magnitude (ascending)
en1=fix(nbl*lowbin*10);                        %number of links to be preserved
en2=fix(en1+nbl);
% en1=round((n^2-n)*lowbin/ud);               %number of links to be preserved
% en2=round((n^2-n)*(lowbin+0.1)/ud);         %number of links to be preserved
E1 = E;
valE1= E1(:,2);
valE1(1:en1)=0;
valE1(en2+1:end)=0;
E1(:,2) = valE1;

for l=1:numel(ind)
    mynode = E1(l,1);
    W(mynode)=E1(l,2);
end

if ud==2                                    %if symmetric matrix
    W=W+W.';                                %reconstruct symmetry
end

end
