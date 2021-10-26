import hail as hl


PLOF_CSQS = ["transcript_ablation", "splice_acceptor_variant",
             "splice_donor_variant", "stop_gained", "frameshift_variant"]

MISSENSE_CSQS = ["stop_lost", "start_lost", "transcript_amplification",
                 "inframe_insertion", "inframe_deletion", "missense_variant"]

SYNONYMOUS_CSQS = ["stop_retained_variant", "synonymous_variant"]

OTHER_CSQS = ["mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"]

# Unused: protein_altering_variant, incomplete_terminal_codon_variant, coding_sequence_variant
# TODO: question, what to do with: "splice_region_variant"
# TODO: question, "missense-damaging" vs "damaging_missense"

def annotate_dbnsfp(mt, vep_path, vep_fields = '/well/lindgren/dpalmer/SAIGE_gene_munging/data/vep_fields.txt'):
    r'''Annotate matrix table with dbNSFP consequence from external VEP file.'''
    print(f'Annotating with VEP file: {vep_path}')
    
    # Open file containing VEP fields
    with open(vep_fields, 'r') as file:
        fields = file.read().strip().split(',')
    ht = hl.import_vcf(vep_path).rename({'info':'vep'}) 
    
    # Add VEP fields by iteration
    for i in range(len(fields)):
        ht = ht.annotate_rows(
            vep=ht.vep.annotate(
                col=ht.vep.CSQ.map(lambda x: (x.split('\\|')[i]))[0]
                ).rename({'col':f'{fields[i]}'})
        )
    
    # Extract various categories annotations and change type
    ht = ht.annotate_rows(vep = ht.vep.annotate(sift_pred = ht.vep.SIFT_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(polyphen2_hdiv_pred = ht.vep.Polyphen2_HDIV_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(polyphen2_hvar_pred = ht.vep.Polyphen2_HVAR_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(cadd_phred_score = hl.parse_float(ht.vep.CADD_phred)))
    ht = ht.annotate_rows(vep = ht.vep.annotate(revel_score = hl.parse_float(ht.vep.REVEL_score)))
    
    # annotate main table
    mt = mt.annotate_rows(dbnsfp = hl.struct())
    mt = mt.annotate_rows(dbnsfp = mt.dbnsfp.annotate(revel_score = ht.index_rows(mt.locus, mt.alleles).vep.revel_score))
    mt = mt.annotate_rows(dbnsfp = mt.dbnsfp.annotate(cadd_phred_score = ht.index_rows(mt.locus, mt.alleles).vep.cadd_phred_score))
    mt = mt.annotate_rows(dbnsfp = mt.dbnsfp.annotate(polyphen2_hdiv_pred = ht.index_rows(mt.locus, mt.alleles).vep.polyphen2_hdiv_pred))
    mt = mt.annotate_rows(dbnsfp = mt.dbnsfp.annotate(polyphen2_hvar_pred = ht.index_rows(mt.locus, mt.alleles).vep.polyphen2_hvar_pred))

    return(mt)


def annotation_case_builder(worst_csq_by_gene_canonical_expr,
                            use_loftee: bool = True,
                            use_polyphen_and_sift: bool = False,
                            use_revel_and_cadd: bool = True):

    case = hl.case(missing_false=True)

    if use_loftee:
        case = (case
                .when(worst_csq_by_gene_canonical_expr.lof == 'HC', 'pLoF')
                .when(worst_csq_by_gene_canonical_expr.lof == 'LC', 'LC'))
    else:
        case = case.when(hl.set(PLOF_CSQS).contains(worst_csq_by_gene_canonical_expr.most_severe_consequence), 'pLoF')
    
    if use_polyphen_and_sift:
        case = (case
                .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence) &
                      (mt.vep.worst_csq_for_variant_canonical.polyphen_prediction == "probably_damaging") &
                      (mt.vep.worst_csq_for_variant_canonical.sift_prediction == "deleterious"), "damaging_missense")
                .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence), "other_missense"))
    else:
        if use_revel_and_cadd:
            case = (case
                .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence) &
                      (mt.vep.worst_csq_for_variant_canonical.polyphen_prediction == "probably_damaging") &
                      (mt.vep.worst_csq_for_variant_canonical.sift_prediction == "deleterious"), "damaging_missense")
                .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence), "other_missense"))
        else:
            case = (case
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_by_gene_canonical_expr.most_severe_consequence) & 
                      (csq_dbnsfp_expr.cadd_phred_score >= 20) & 
                      (csq_dbnsfp_expr.revel_score >= 0.6), "damaging_missense")
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_by_gene_canonical_expr.most_severe_consequence), "other_missense"))

    case = case.when(hl.set(SYNONYMOUS_CSQS).contains(worst_csq_by_gene_canonical_expr.most_severe_consequence), 'synonymous')
    case = case.when(hl.set(OTHER_CSQS).contains(worst_csq_by_gene_canonical_expr.most_severe_consequence), 'non-coding')
    
    return case.or_missing()


def create_gene_map_ht(ht, check_gene_contigs=False):
    from gnomad.utils.vep import process_consequences

    ht = process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.annotate(
        variant_id=ht.locus.contig + ':' + hl.str(ht.locus.position) + '_' + ht.alleles[0] + '/' + ht.alleles[1],
        annotation=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical))
    if check_gene_contigs:
        gene_contigs = ht.group_by(
            gene_id=ht.vep.worst_csq_by_gene_canonical.gene_id,
            gene_symbol=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
        ).aggregate(
            contigs=hl.agg.collect_as_set(ht.locus.contig)
        )
        assert gene_contigs.all(hl.len(gene_contigs.contigs) == 1)

    gene_map_ht = ht.group_by(
        gene_id=ht.vep.worst_csq_by_gene_canonical.gene_id,
        gene_symbol=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
    ).partition_hint(100).aggregate(
        interval=hl.interval(
            start=hl.locus(hl.agg.take(ht.locus.contig, 1)[0], hl.agg.min(ht.locus.position)),
            end=hl.locus(hl.agg.take(ht.locus.contig, 1)[0], hl.agg.max(ht.locus.position))
        ),
        variants=hl.agg.group_by(ht.annotation, hl.agg.collect(ht.variant_id)),
    )
    return gene_map_ht


def post_process_gene_map_ht(gene_ht):
    groups = ['pLoF', 'damaging_missense|LC', 'pLoF|damaging_missense|LC', 'pLoF|damaging_missense',  'damaging_missense', 'other_missense', 'synonymous']
    variant_groups = hl.map(lambda group: group.split('\\|').flatmap(lambda csq: gene_ht.variants.get(csq)), groups)
    gene_ht = gene_ht.transmute(
        variant_groups=hl.zip(groups, variant_groups)
    ).explode('variant_groups')
    gene_ht = gene_ht.transmute(annotation=gene_ht.variant_groups[0], variants=hl.sorted(gene_ht.variant_groups[1]))
    gene_ht = gene_ht.key_by(start=gene_ht.interval.start)
    return gene_ht.filter(hl.len(gene_ht.variants) > 0)
