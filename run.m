topdir = '/path/to/ApoHP';
srcdir = [topdir '/src/'];
addpath(srcdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% survey_hairpins
% (takes ~1 day to run in loop mode, or ~1 hr to run in parallelized mode)

refdir = [topdir '/ref/hg19'];
datdir = [topdir '/data/processed/'];

ApoHP('survey_hairpins',refdir,datdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% analyze_mutations
% (takes ~4 hrs to run)

refdir = [topdir '/ref/hg19'];
initdir = [topdir '/data/init/'];
datdir = [topdir '/data/processed/'];
figdir = [topdir '/figs'];

allCsfile = [datdir '/allCs.mat'];
ttypefile = [initdir '/tumortype_longnames_and_colors.txt'];
mutationfile = [initdir '/mutations.wgs.txt'];
sigfile = [initdir '/sigProfiler_SBS_signatures.WGS.v3.0.mat'];
k = 8;
wxs_mutationfile = [initdir '/mutations.wxs.txt'];
wxs_sigfile = [initdir '/sigProfiler_SBS_signatures.WXS.v3.0.mat'];
wxs_k = 12;

ApoHP('analyze_mutations',refdir, mutationfile, ttypefile, sigfile, k, allCsfile, datdir, figdir, wxs_mutationfile, wxs_sigfile, wxs_k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% superimpose_data
% (takes ~30 min to run)

refdir = [topdir '/ref/hg19'];
initdir = [topdir '/data/init/'];
datdir = [topdir '/data/processed/'];
figdir = [topdir '/figs'];

allCsfile = [datdir '/allCs.mat'];
fig3file = [datdir '/figure3.mat'];
ttypefile = [initdir '/tumortype_longnames_and_colors.txt'];
mutationfile = [initdir '/mutations.wxs.500mut.txt'];
sigfile = [initdir '/sigProfiler_SBS_signatures.WXS.v3.0.mat'];
k = 12;

ApoHP('superimpose_data',refdir, mutationfile, ttypefile, sigfile, k, allCsfile, fig3file, datdir, figdir);




