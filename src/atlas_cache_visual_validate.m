function atlas_cache_visual_validate(S, cfg, outdir, cacheInfo)
%ATLAS_CACHE_VISUAL_VALIDATE  Step 5 visuals: cache stats summary.

if ~cfg.io.save_figs
    return;
end

st = atlas_cache_stats(S, cfg);

f = figure('Visible', cfg.io.fig_visible); clf;
vals = [st.rows_FRS_upper, st.rows_BRS_upper, st.rows_FRS_lower, st.rows_BRS_lower];
bar(vals);
set(gca,'XTickLabel',{'FRS_u','BRS_u','FRS_l','BRS_l'});
ylabel('rows');
title(sprintf('Step5 cache stats: %s | hit=%d | %s', S.name, cacheInfo.hit, cacheInfo.hash));
grid on;

drawnow;
safeName = sanitize_fname(S.name);
saveas(f, fullfile(outdir, sprintf('step5_%s_cache_stats.png', safeName)));
close(f);
end
