function ApoHP(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% APOBEC hairpins analysis tool
% (C) 2021 Adam Langenbucher and Michael S. Lawrence
% Nature Communications https://www.nature.com/articles/s41467-021-21891-0
%
% ApoHP <command> <args>
%
% command = { survey_hairpins, analyze_mutations, superimpose_data }
%
% apoHP survey_hairpins <refdir> <outdir> [<blockno>]
% --> reads genomic reference data from <refdir>
% --> scans each genomic block for hairpins (saves temporary block files in outdir>
% --> collects all blocks
% --> NOTE: <blockno> can be specified if running in parallel (specify -1 to perform "gather" step after all blocks have finished)
% --> writes <outdir>/allCs.mat                with list of all Cs in genome and their hairpin characteristics
% --> writes <outdir>/best_hairpin_sites.mat   with list of best hairpin Cs (3/3 ss>=10, 4/[4 or 5] ss>=12)
% --> writes <outdir>/best_hairpin_sites.bed   (same but in text format)
%
% ApoHP analyze_mutations <refdir> <mutationfile> <ttypefile> <sigfile> <k> <allCsfile> <datoutdir> <figoutdir> <wxs_mutationfile> <wxs_sigfile> <wxs_k>
% --> reads mutations from input file
% --> runs NMF to discover <k> mutation signatures
% --> maps mutations to allCs file
% --> identifies 10%, 50%, and 90% APOBEC cohorts, including or excluding MSUPE patients
% --> tabulates number of mutated patients at each C site
% --> calculates YTCA/RTCA and TpC/VpC %hairpin statistics
% --> writes text files in <datoutdir> of:
%        1. patients and their signature loadings and APOBEC metrics
%        2. data for figures
% --> writes figure PDFs in <figoutdir>
% --> writes <datoutdir>/full_data.mat                       with all data
% --> writes <datoutdir>/full_data.no_wgs_mutations.mat      with all data except X.mut
%
% ApoHP superimpose_data <refdir> <mutationfile> <ttypefile> <sigfile> <k> <allCsfile> <fig3file> <datoutdir> <figoutdir>
% --> loads user-provided mutation data from <newmutfile>
% --> combines user-provided mutation set with pan-cancer mutation set
% --> runs NMF to determine frac_APOBEC mutational signature loading
% --> maps user-provided mutation set to annotated C:G basepair object
% --> computes APOBEC3A/3B character for each sample in user-defined cohort
% --> computes TpC/VpC hairpin mutation prevalence for each sample in user-defined cohort
% --> writes figures in <figoutdir>
% --> writes output data tables in <datoutdir>

fprintf('ApoHP:  APOBEC hairpins analysis tool\n');  
fprintf('(C) 2020 Adam Langenbucher and Michael S. Lawrence\n\n');

if nargin<1, error('usage: ApoHP <command> <args>'); end
command = varargin{1};

fprintf('command: %s\n\n',command)

if strcmp(command, 'survey_hairpins')

  if nargin<3, error('requires two arguments: refdir and outdir'); end
  refdir = varargin{2};
  outdir = varargin{3};
  blockno = nan;
  if nargin>=4
    blockno = num2str(varargin{4});
    if ~(blockno>=1 || blockno==-1 && blockno==round(blockno)), error('blockno not valid (should be {1,2,...,#blocks} or -1 for gather'); end
  end 
  if nargin>4, error('Too many input arguments'); end

  X = load_genome_info(refdir);
  fprintf('Number of 10Mb blocks in genome: %d\n',slength(X.block));
  
  if isnan(blockno)                                % non-parallel mode: run as loop
    fprintf('Surveying genome for hairpins (%d blocks).\n',slength(X.block));
    for blockno=1:slength(X.block)
      survey_hairpins(refdir,outdir,blockno);
    end
    survey_hairpins_gather(refdir,outdir);
  elseif blockno>=1 && blockno<=slength(X.block)   % parallel mode: "scatter"
    survey_hairpins(refdir,outdir,blockno);
  elseif blockno==-1                               % parallel mode: "gather"
    survey_hairpins_gather(refdir,outdir);
  else
    error('blockno %d out of range',blockno);
  end
  
elseif strcmp(command, 'analyze_mutations')

  args = varargin(2:end);
  analyze_mutations(args{:});

elseif strcmp(command, 'superimpose_data')

  args = varargin(2:end);
  superimpose_data(args{:});

else

  error('Unknown command.');

end

fprintf('\nApoHP finished successfully.\n');
