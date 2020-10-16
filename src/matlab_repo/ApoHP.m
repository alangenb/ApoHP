function ApoHP(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% APOBEC hairpins tool
%
% ApoHP <command> <args>
%
% command = { survey_hairpins, analyze_mutations }
%
% apoHP survey_hairpins <refdir> <outdir>
% --> reads genomic reference data from <refdir>
% --> scans each genomic block for hairpins (saves temporary block files in outdir>
% --> collects all blocks
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
% --> writes PDFs in <figoutdir> of:
%        1. lego plots of mutation signatures
%        2. figures from paper
% --> writes <datoutdir>/full_data.mat                       with all data
% --> writes <datoutdir>/full_data.no_mutation_list.mat      with all data except X.mut

command = varargin{1};

if strcmpi(command, 'survey_hairpins')

  refdir = varargin{2};
  outdir = varargin{3};
  if nargin>3, error('Too many input arguments'); end

  X = load_genome_info(refdir);
  fprintf('Surveying genome for hairpins (%d blocks).\n',slength(X.block));
  for blockno=1:slength(X.block)
    survey_hairpins(refdir,outdir,blockno);
  end
  survey_hairpins_gather(refdir,outdir);

elseif strcmpi(command, 'analyze_mutations')

  args = varargin(2:end);
  analyze_mutations(args{:});

else

  error('Unknown command "%s"', command);

end

