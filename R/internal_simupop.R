main_12_simupop_api <- function(args = commandArgs(trailingOnly = TRUE)) {

# R API wrapper around SimuPOP via reticulate (inline Python)

args <- args
usage <- function() {
  cat("Usage: 12_simupop_api.R --config <path> --out_prefix <path>\n")
  quit(status = 1)
}

get_arg <- function(flag, default = NULL) {
  if (!(flag %in% args)) return(default)
  idx <- match(flag, args)
  if (idx == length(args)) return(default)
  args[idx + 1]
}

cfg_path <- get_arg("--config")
out_prefix <- get_arg("--out_prefix")
if (is.null(cfg_path) || is.null(out_prefix)) usage()

if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required for SimuPOP configuration.")
}
if (!requireNamespace("reticulate", quietly = TRUE)) {
  stop("Package 'reticulate' is required for SimuPOP integration.")
}

cfg <- jsonlite::fromJSON(cfg_path)

python_hook <- cfg$python_hook
if (!is.null(python_hook) && nzchar(python_hook)) {
  reticulate::py_run_string(python_hook)
}

if ("simitall" %in% reticulate::conda_list()$name) {
  reticulate::use_condaenv("simitall", required = TRUE)
}

py_code <- "
import random
import simuPOP as sim


def _make_parent_chooser(cfg):
    chooser = (cfg.get('chooser') or 'RandomParentsChooser').strip()
    mapping = {
        'SequentialParentChooser': sim.SequentialParentChooser,
        'SequentialParentsChooser': sim.SequentialParentsChooser,
        'RandomParentChooser': sim.RandomParentChooser,
        'RandomParentsChooser': sim.RandomParentsChooser,
        'PolyParentsChooser': sim.PolyParentsChooser,
        'CombinedParentsChooser': sim.CombinedParentsChooser,
        'PyParentsChooser': sim.PyParentsChooser,
    }
    if chooser not in mapping:
        raise ValueError(f'Unknown parent chooser: {chooser}')
    if chooser == 'PyParentsChooser':
        gen_name = cfg.get('generator_name')
        if not gen_name or gen_name not in globals():
            raise ValueError('PyParentsChooser requires generator_name from python_hook globals')
        return mapping[chooser](globals()[gen_name])
    if chooser == 'CombinedParentsChooser':
        fc = cfg.get('fatherChooser') or {}
        mc = cfg.get('motherChooser') or {}
        father = _make_parent_chooser(fc)
        mother = _make_parent_chooser(mc)
        return mapping[chooser](father, mother, allowSelfing=cfg.get('allowSelfing', True))
    kwargs = {}
    for k in ['replacement', 'selectionField', 'sexChoice', 'polySex', 'polyNum']:
        if k in cfg and cfg[k] is not None:
            kwargs[k] = cfg[k]
    return mapping[chooser](**kwargs) if kwargs else mapping[chooser]()


def _make_offspring_generator(cfg):
    gen = (cfg.get('generator') or 'OffspringGenerator').strip()
    if gen == 'ControlledOffspringGenerator':
        freq_name = cfg.get('freqFunc_name')
        if not freq_name or freq_name not in globals():
            raise ValueError('ControlledOffspringGenerator requires freqFunc_name from python_hook globals')
        loci = cfg.get('loci', [])
        alleles = cfg.get('alleles', [])
        return sim.ControlledOffspringGenerator(loci=loci, alleles=alleles, freqFunc=globals()[freq_name])
    if gen != 'OffspringGenerator':
        raise ValueError(f'Unknown offspring generator: {gen}')
    num = cfg.get('numOffspring', 1)
    sex_mode = cfg.get('sexMode')
    kwargs = {'numOffspring': num}
    if sex_mode is not None:
        kwargs['sexMode'] = sex_mode
    return sim.OffspringGenerator(ops=[sim.MendelianGenoTransmitter()], **kwargs)


def _make_mating_scheme(cfg):
    scheme = (cfg.get('scheme') or 'RandomMating').strip()
    size = cfg.get('offspring') if cfg.get('offspring') is not None else cfg.get('subPopSize')

    # Directly supported concrete schemes
    mapping = {
        'MatingScheme': sim.MatingScheme,
        'CloneMating': sim.CloneMating,
        'RandomSelection': sim.RandomSelection,
        'RandomMating': sim.RandomMating,
        'MonogamousMating': sim.MonogamousMating,
        'PolygamousMating': sim.PolygamousMating,
        'HaplodiploidMating': sim.HaplodiploidMating,
        'SelfMating': sim.SelfMating,
        'HermaphroditicMating': sim.HermaphroditicMating,
        'ControlledRandomMating': sim.ControlledRandomMating,
    }

    if scheme in mapping:
        kwargs = {}
        if size is not None:
            kwargs['subPopSize'] = size
        for k in ['numOffspring', 'sexMode', 'selectionField', 'polySex', 'polyNum', 'replacement', 'allowSelfing', 'ops']:
            if k in cfg and cfg[k] is not None:
                kwargs[k] = cfg[k]
        return mapping[scheme](**kwargs) if kwargs else mapping[scheme]()

    if scheme == 'HomoMating':
        chooser = _make_parent_chooser(cfg)
        generator = _make_offspring_generator(cfg)
        return sim.HomoMating(chooser, generator,
                              subPopSize=size,
                              subPops=cfg.get('subPops', sim.ALL_AVAIL),
                              weight=cfg.get('weight', 0))

    if scheme == 'HeteroMating':
        schemes_cfg = cfg.get('matingSchemes') or []
        if not schemes_cfg:
            raise ValueError('HeteroMating requires matingSchemes list')
        schemes = [_make_mating_scheme(sc) for sc in schemes_cfg]
        return sim.HeteroMating(schemes,
                                subPopSize=size,
                                shuffleOffspring=cfg.get('shuffleOffspring', True),
                                weightBy=cfg.get('weightBy', sim.ANY_SEX))

    if scheme == 'ConditionalMating':
        cond = cfg.get('cond', True)
        if_cfg = cfg.get('ifMatingScheme')
        else_cfg = cfg.get('elseMatingScheme')
        if if_cfg is None or else_cfg is None:
            raise ValueError('ConditionalMating requires ifMatingScheme and elseMatingScheme')
        return sim.ConditionalMating(cond, _make_mating_scheme(if_cfg), _make_mating_scheme(else_cfg))

    if scheme == 'PedigreeMating':
        ped_name = cfg.get('pedigree_name')
        ops_name = cfg.get('ops_name')
        if not ped_name or ped_name not in globals():
            raise ValueError('PedigreeMating requires pedigree_name from python_hook globals')
        ops = globals()[ops_name] if ops_name and ops_name in globals() else []
        return sim.PedigreeMating(globals()[ped_name], ops)

    raise ValueError(f'Unknown mating scheme: {scheme}')

def _read_fasta_multi(path):
    ids, seqs = [], []
    with open(path) as fh:
        cur_id = None
        cur_seq = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if cur_id is not None:
                    ids.append(cur_id)
                    seqs.append(''.join(cur_seq).upper())
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            ids.append(cur_id)
            seqs.append(''.join(cur_seq).upper())
    return ids, seqs


def _panel_to_loci(seqs, max_loci=None):
    L = len(seqs[0])
    for s in seqs:
        if len(s) != L:
            raise ValueError('All haplotypes must be same length')
    var_positions = []
    alleles_per_pos = []
    for i in range(L):
        alleles = sorted(set(s[i] for s in seqs))
        if len(alleles) > 1:
            var_positions.append(i + 1)
            alleles_per_pos.append(alleles)
    if max_loci and len(var_positions) > max_loci:
        idxs = [int(i * len(var_positions) / max_loci) for i in range(max_loci)]
        var_positions = [var_positions[i] for i in idxs]
        alleles_per_pos = [alleles_per_pos[i] for i in idxs]
    return var_positions, alleles_per_pos


def _hap_to_geno(hap_seq, var_positions, allele_maps):
    out = []
    for pos, amap in zip(var_positions, allele_maps):
        allele = hap_seq[pos - 1]
        out.append(amap.get(allele, 0))
    return out


def _export_genotypes(pop, out_prefix):
    geno_path = f'{out_prefix}.genotypes.tsv'
    meta_path = f'{out_prefix}.meta.tsv'
    loci = pop.totNumLoci()
    ploidy = pop.ploidy()
    with open(geno_path, 'w') as gf:
        gf.write('sample\t' + '\t'.join([f'locus_{i+1}' for i in range(loci)]) + '\n')
        for ind in pop.individuals():
            sid = ind.info('ind_id') if 'ind_id' in pop.infoFields() else 'NA'
            row = []
            for l in range(loci):
                alleles = [str(ind.allele(l, ploidy=p)) for p in range(ploidy)]
                row.append('|'.join(alleles))
            gf.write(str(sid) + '\t' + '\t'.join(row) + '\n')
    with open(meta_path, 'w') as mf:
        mf.write('sample\n')
        for ind in pop.individuals():
            sid = ind.info('ind_id') if 'ind_id' in pop.infoFields() else 'NA'
            mf.write(str(sid) + '\n')


def _export_vcf(pop, out_prefix):
    vcf_path = f'{out_prefix}.vcf'
    loci = pop.totNumLoci()
    ploidy = pop.ploidy()
    samples = [str(ind.info('ind_id')) for ind in pop.individuals()]
    with open(vcf_path, 'w') as vf:
        vf.write('##fileformat=VCFv4.2\n')
        vf.write(f'##contig=<ID=chr1,length={loci}>\n')
        vf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n')
        vf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')
        for l in range(loci):
            row = ['chr1', str(l+1), f'var{l+1}', '0', '1', '.', 'PASS', '.', 'GT']
            gts = []
            for ind in pop.individuals():
                alleles = [str(ind.allele(l, ploidy=p)) for p in range(ploidy)]
                gt = alleles[0] if ploidy == 1 else '/'.join(alleles[:2])
                gts.append(gt)
            row.extend(gts)
            vf.write('\t'.join(row) + '\n')


def run_simupop(config, out_prefix):
    pop_cfg = config.get('population', {})
    mating_cfg = config.get('mating', {})
    init_cfg = config.get('init', {})
    preset = config.get('preset')
    gens = int(config.get('generations', 1))

    size = pop_cfg.get('size', 100)
    ploidy = pop_cfg.get('ploidy', 2)
    loci = pop_cfg.get('loci', [100])
    chrom_names = pop_cfg.get('chromNames', [])
    loci_pos = pop_cfg.get('lociPos', [])
    loci_names = pop_cfg.get('lociNames', [])
    allele_names = pop_cfg.get('alleleNames', [])
    info_fields = pop_cfg.get('infoFields', ['ind_id'])

    init_from_fasta = init_cfg.get('from_fasta')
    max_loci = init_cfg.get('max_loci')
    sample_haplotypes = init_cfg.get('sample_haplotypes', True)
    track_hap_ids = init_cfg.get('track_hap_ids', False)

    if init_from_fasta:
        hap_ids, hap_seqs = _read_fasta_multi(init_from_fasta)
        var_positions, alleles_per_pos = _panel_to_loci(hap_seqs, max_loci=max_loci)
        if not var_positions:
            raise ValueError('No variable sites found in haplotype panel')
        loci = [len(var_positions)]
        loci_pos = var_positions
        allele_names = alleles_per_pos
        chrom_names = ['chr1']

    pop = sim.Population(size=size, ploidy=ploidy, loci=loci, chromNames=chrom_names,
                         lociPos=loci_pos, lociNames=loci_names, alleleNames=allele_names,
                         infoFields=info_fields)

    for i, ind in enumerate(pop.individuals()):
        ind.setInfo(i+1, 'ind_id')

    if preset == 'F1':
        mating_cfg = {'scheme': 'RandomMating', 'offspring': pop.popSize()}
        gens = 1
    elif preset == 'F2':
        mating_cfg = {'scheme': 'SelfMating', 'offspring': pop.popSize()}
        gens = 2
    elif preset == 'MAGIC':
        mating_cfg = {'scheme': 'RandomMating', 'offspring': pop.popSize()}
        gens = 5
    elif preset == 'NAM':
        mating_cfg = {'scheme': 'RandomMating', 'offspring': pop.popSize()}
        gens = 3

    if init_from_fasta:
        if track_hap_ids:
            pop.addInfoFields(['hap1', 'hap2'], init=0)
        allele_maps = [{a:i for i,a in enumerate(alleles)} for alleles in allele_names]
        for ind in pop.individuals():
            hap_idxs = [random.randrange(len(hap_seqs)) for _ in range(ploidy)] if sample_haplotypes else [0 for _ in range(ploidy)]
            geno = []
            for hi in hap_idxs:
                geno.extend(_hap_to_geno(hap_seqs[hi], loci_pos, allele_maps))
            ind.setGenotype(geno)
            if track_hap_ids and ploidy >= 2:
                ind.setInfo(hap_ids[hap_idxs[0]], 'hap1')
                ind.setInfo(hap_ids[hap_idxs[1]], 'hap2')

    mating_scheme = _make_mating_scheme(mating_cfg)
    sim.Simulator(pop).evolve(matingScheme=mating_scheme, gen=gens)

    _export_genotypes(pop, out_prefix)
    if config.get('export_vcf', False):
        _export_vcf(pop, out_prefix)

    return {'genotypes': f'{out_prefix}.genotypes.tsv', 'meta': f'{out_prefix}.meta.tsv'}
"

reticulate::py_run_string(py_code)

res <- reticulate::py$run_simupop(cfg, out_prefix)
cat("Wrote:\n")
cat("  ", res$genotypes, "\n", sep = "")
cat("  ", res$meta, "\n", sep = "")
}
