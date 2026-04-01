function charDrawGaps(ax, gapIn, Charcoal)
% charDrawGaps   Draw grey patch overlays at missing-value gaps.
%
%   charDrawGaps(ax, gapIn, Charcoal)
%
%   Draws a grey filled patch over each gap interval on axes ax.
%   gapIn is an [nGaps x 2] matrix of indices into Charcoal.ybp.
%   If ax is empty, uses the current axes.

if isempty(ax)
    ax = gca;
end
yl = get(ax, 'ylim');
for gi = 1:size(gapIn, 1)
    xg = [Charcoal.ybp(gapIn(gi,1))  Charcoal.ybp(gapIn(gi,1)) ...
          Charcoal.ybp(gapIn(gi,2))  Charcoal.ybp(gapIn(gi,2))];
    patch(xg, [yl(1) yl(2) yl(2) yl(1)], 'w', ...
        'edgecolor', [.75 .75 .75], 'facecolor', [.75 .75 .75])
end

end
