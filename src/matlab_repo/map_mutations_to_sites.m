function site_idx = map_mutations_to_sites(mut,site)
% site_idx = map_mutations_to_sites(mut,site)
% --> returns site_idx, an index into site list

demand_fields(mut,{'chr','pos'});
demand_fields(site,{'chr','pos'});
  
site_idx = ismembc2(1e9*double(mut.chr)+double(mut.pos),1e9*double(site.chr)+double(site.pos));

if all(isnan(site_idx)), fprintf('Warning: Nothing mapped.\n'); end


