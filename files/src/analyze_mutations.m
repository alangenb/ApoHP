function X = analyze_mutations(refdir, mutationfile, ttypefile, sigfile, k, allCsfile, datoutdir, figoutdir, wxs_mutationfile, wxs_sigfile, wxs_k)
% analyze_mutations(refdir, mutationfile, ttypefile, sigfile, k, allCsfile, datoutdir, figoutdir, [wxs_mutationfile, wxs_sigfile, wxs_k])
% --> reads mutations from input file
% --> runs NMF to discover <k> mutation signatures
% --> maps mutations to allCs file
% --> identifies 10%, 50%, and 90% APOBEC cohorts, including or excluding MSUPE patients
% --> tabulates number of mutated patients at each C site
% --> calculates YTCA/RTCA and TpC/VpC %hairpin statistics
% --> writes allCs.<mutation_file>.mat with all analyses (but no X.mut)
% --> writes PDFs:
%        1. lego plots of mutation signatures
%        2. Figure 3-type plots
% --> writes text files of:
%        1. patients and their signature loadings and APOBEC metrics

if nargin<8, error('requires: refdir, mutationfile, ttypefile, sigfile, k, allCsfile, datoutdir, figoutdir');
elseif nargin==8 || nargin==11
  demand_files({refdir,mutationfile,ttypefile,sigfile,allCsfile});
  if ~isnumeric(k), k = str2double(k); end
  if ~(k>=2 & k<=1000), error('k should be 2-1000'); end
  if nargin==11
    demand_files({wxs_mutationfile,wxs_sigfile});
    if ~isnumeric(wxs_k), wxs_k = str2double(wxs_k); end
    if ~(wxs_k>=2 & wxs_k<=1000), error('wxs_k should be 2-1000'); end
    wxs_provided = true;
  else
    wxs_provided = false;
  end
else
  error('inputs should be: refdir, mutationfile, ttypefile, sigfile, k, allCsfile, datoutdir, figoutdir, [wxs_mutationfile, wxs_sigfile, wxs_k]');
end

% LOAD MUTATION DATA
% will save as compact MAT file for quicker re-loading
ede(datoutdir);
mutname = regexprep(mutationfile,'^.*\/([^/]+)$','$1');
mutationmat = [datoutdir '/' mutname '.mat'];
X = load_mutation_data(mutationfile,refdir,mutationmat,ttypefile);

% SIGNATURE ANALYSIS
randseed=1234;
X = analyze_mutation_signatures(X,k,randseed,sigfile);

% LOAD allCs file
tmp=load(allCsfile,'X');tmp=tmp.X;
fs = fieldnames(tmp);
for i=1:length(fs),f=fs{i};X.(f)=tmp.(f);end
clear tmp

% MAP MUTATIONS TO Cs
X.mut.site_idx = map_mutations_to_sites(X.mut,X.site);

% COMPUTE %HP METRICS
X = compute_pcthp_metrics(X);

% IDENTIFY RECURRENT APOBEC MUTATIONS
p = (X.pat.frac_apobec>=0.1);
X.site.wgs_apo10_ct = histc(X.mut.site_idx(p(X.mut.pat_idx)),1:slength(X.site));
X.site.wgs_apo10_pct = 100*X.site.wgs_apo10_ct/sum(p);

ci=4;
X.site.wgs_apo10_ct = histc(X.mut.site_idx(p(X.mut.pat_idx)),1:slength(X.site));
X.site.wgs_apo10_pct = 100*X.site.wgs_apo10_ct/sum(p);

% PROCESS WXS DATA if provided
if wxs_provided
  fprintf('Processing WXS data.\n');
  wxs_mutname = regexprep(wxs_mutationfile,'^.*\/([^/]+)$','$1');
  wxs_mutationmat = [datoutdir '/' wxs_mutname '.mat'];
  X.wxs = load_mutation_data(wxs_mutationfile,refdir,wxs_mutationmat,ttypefile);
  randseed=1234;
  X.wxs = analyze_mutation_signatures(X.wxs,wxs_k,randseed,wxs_sigfile);
  X.wxs.mut.site_idx = map_mutations_to_sites(X.wxs.mut,X.site);
  p_wxs = (X.wxs.pat.frac_apobec>=0.1);
  X.site.wxs_apo10_ct = histc(X.wxs.mut.site_idx(p_wxs(X.wxs.mut.pat_idx)),1:slength(X.site));
  X.site.wxs_apo10_pct = 100*X.site.wxs_apo10_ct/sum(p_wxs);
end

% SAVE PATIENTS LIST
pat = X.pat;
outname = [datoutdir '/patients.mat'];
save(outname,'pat','-v7.3');

% TABULATE CLASSES OF Cs and their mutation counts
X = tabulate_sites(X);

% QUANTITATIVE MODELING
ci=1;
X = apobec_quantitative_modeling(X,ci);

% SAVE FULL OBJECT
outname = [datoutdir '/full_data.mat'];
save(outname,'X','-v7.3')

% SAVE FULL OBJECT without WGS mutation list (can distribute this publicly without violating TCGA WGS rules)
X = rmfield(X,'mut');
outname = [datoutdir '/full_data.no_wgs_mutations.mat'];
save(outname,'X','-v7.3')

% MAKE DATA FILES FOR FIGURES
make_data_files_for_figures(X,datoutdir);

% MAKE FIGURES
make_figures(datoutdir,figoutdir);


