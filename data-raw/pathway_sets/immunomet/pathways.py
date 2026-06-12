#!/usr/bin/env python3
"""
Immunometabolism set membership rules — O'Neill framework.

Based primarily on:
  O'Neill LAJ, Kishton RJ, Rathmell J. (2016). A guide to immunometabolism for
  immunologists. Nat Rev Immunol 16(9):553-565. doi:10.1038/nri.2016.70

Supporting references (for context on specific metabolites):
  Pearce EL, Pearce EJ. (2013). Metabolic pathways in immune cell activation
  and quiescence. Immunity 38(4):633-643.
  
  Buck MD, O'Sullivan D, Pearce EL. (2015). T cell metabolism drives immunity.
  J Exp Med 212(9):1345-1360.
  
  Mogilenko DA, Sergushichev A, Artyomov MN. (2023). Systems Immunology
  Approaches to Metabolism. Annu Rev Immunol 41:317-342.

The original O'Neill 2016 framework enumerates six core pathways. The pentose
phosphate pathway is excluded because the Biocrates Quant 1000 panel does not
cover its canonical metabolites (G6P, 6-P-gluconate, ribose-5-P, NADPH).
The remaining five pathways are defined below.
"""

# =============================================================================
# Canonical metabolites per pathway
# =============================================================================
#
# Each pathway lists the metabolites that O'Neill 2016 (and supporting refs)
# treat as core pathway members. Names here match the Biocrates Quant 1000
# `shortname` or `fullname` field where possible. Metabolites that are
# canonical to the pathway but ABSENT from the Biocrates panel are listed
# in comments so that coverage gaps are documented.

# -----------------------------------------------------------------------------
# GLYCOLYSIS
# -----------------------------------------------------------------------------
# O'Neill 2016 Figure 1: glucose → pyruvate → lactate, with HK, PFK, LDH as
# key enzymes. In activated immune cells (M1 macrophages, effector T cells):
# aerobic glycolysis (Warburg-like) with enhanced lactate output.
#
# Biocrates Q1000 coverage: Glucose, Pyruvic acid (in Carboxylic acids),
# Lac (Lactic acid, in Carboxylic acids). Upstream intermediates (G6P,
# F1,6BP, PEP) are NOT in the panel.
GLYCOLYSIS_METABOLITES = {
    'Glucose',       # panel class: Sugars
    'Pyruvic acid',  # panel class: Carboxylic acids
    'Lac',           # panel class: Carboxylic acids (Lactic acid)
}
# Absent from panel (acknowledged gaps):
#   - G6P, F6P, F1,6BP, GAP, DHAP, 1,3-BPG, 3-PG, 2-PG, PEP (no panel coverage)

# -----------------------------------------------------------------------------
# TCA CYCLE
# -----------------------------------------------------------------------------
# O'Neill 2016: TCA intermediates with "broken" TCA in M1 macrophages
# (citrate & succinate accumulation, itaconate production via ACOD1/IRG1).
# Succinate drives HIF-1α stabilization; itaconate is anti-inflammatory via
# SDH inhibition and Nrf2 activation.
#
# Biocrates Q1000: Citric acid, Isocitric acid, AconAcid (Aconitate),
# a-Ketoglutaric acid, Suc (Succinic acid), Fumaric acid, Malic acid,
# Oxaloacetic acid, 2-OH-Glutaric acid, Itaconic acid.
TCA_METABOLITES = {
    'Citric acid',
    'Isocitric acid',
    'AconAcid',           # Aconitate / cis-aconitate
    'a-Ketoglutaric acid',
    'Suc',                # Succinic acid
    'Fumaric acid',
    'Malic acid',
    'Oxaloacetic acid',
    'Itaconic acid',      # ACOD1 product; immunometabolism core
    '2-OH-Glutaric acid', # oncometabolite; IDH mutation context
}

# -----------------------------------------------------------------------------
# FATTY ACID OXIDATION (FAO)
# -----------------------------------------------------------------------------
# O'Neill 2016: FAO preferred by memory T cells, M2 macrophages, Tregs.
# Carnitine shuttle (CPT1/CPT2) moves acyl groups into mitochondria.
# Acylcarnitines are the diagnostic signature of FAO flux.
# Dicarboxylic acids arise from peroxisomal ω-oxidation and are indicators
# of incomplete / alternative fatty acid oxidation.
#
# Biocrates Q1000: Full acylcarnitine class (40 metabolites) + Free carnitine
# + Dicarboxylic acids class.
FAO_METABOLITE_CLASSES = {
    'Acylcarnitines',       # includes C0 (free carnitine), C2-C26 acylcarnitines
    'Dicarboxylic acids',   # products of ω-oxidation / incomplete β-oxidation
}
# Individual non-class metabolites that belong here:
FAO_INDIVIDUAL = set()  # none additional; all FAO-relevant in above classes

# Explicit exclusion: the Biocrates "Dicarboxylic acids" class chemically contains
# several TCA-cycle intermediates (succinate, fumarate, malate, oxaloacetate,
# α-ketoglutarate, 2-hydroxyglutarate, itaconate). These are TCA-specific and
# are NOT products of fatty acid β- or ω-oxidation. We therefore exclude them
# from FAO so that Dicarboxylic acids matched in FAO are only the true
# ω-oxidation / alternative-FA-oxidation products (adipic, suberic, DiCA 12:0, etc.)
# and BCAA-catabolism byproducts that functionally pair with FAO.
FAO_EXCLUDE = {
    '2-OH-Glutaric acid',
    'a-Ketoglutaric acid',
    'Fumaric acid',
    'Itaconic acid',
    'Malic acid',
    'Oxaloacetic acid',
    'Suc',
    # (Citric acid, Isocitric acid, AconAcid are in Tricarboxylic acids class,
    #  which is not in FAO_METABOLITE_CLASSES, so they're already excluded.)
}

# -----------------------------------------------------------------------------
# FATTY ACID SYNTHESIS (FAS)
# -----------------------------------------------------------------------------
# O'Neill 2016: De novo FAS required for activated DC membrane expansion,
# Th17 differentiation (ACC1 dependency), effector T cell function.
# ACC1 → malonyl-CoA → palmitate; SCD for desaturation; then esterified into
# complex lipids. In plasma metabolomics, de novo lipogenesis is inferred
# from saturated fatty acids, monounsaturated fatty acids (SCD index), and
# newly synthesized glycerolipids/phospholipids.
#
# Biocrates Q1000: Fatty acids class (individual FAs), Triglycerides,
# Diglycerides, Monoglycerides, Cholesteryl esters,
# Phosphatidylcholines, Lysophosphatidylcholines (choline-containing
# phospholipids are the primary readout of de novo lipogenesis in plasma).
FAS_METABOLITE_CLASSES = {
    'Fatty acids',                   # individual free fatty acids
    'Triglycerides',                 # TG storage / assembly
    'Diglycerides',                  # DG intermediates
    'Monoglycerides',                # MG intermediates
    'Cholesteryl esters',            # cholesterol esterification products
    'Phosphatidylcholines',          # PC — primary de novo product
    'Lysophosphatidylcholines',      # LPC — derived from PC, informs PC pool
}
FAS_INDIVIDUAL = set()

# Note on overlap with FAO: Fatty acids class appears in both FAS (as the
# product of synthesis) and FAO (as the substrate for oxidation). In plasma,
# free fatty acids reflect net flux of both processes and cannot be
# cleanly assigned to one. O'Neill 2016 treats FAS and FAO as distinct
# pathways; we honor that by allowing Fatty acids to be in both sets.
# However, to avoid double-counting the substrate-vs-product ambiguity,
# only acylcarnitines + dicarboxylic acids (FAO-specific products) go in FAO,
# while Fatty acids + glycerolipids + choline-lipids go in FAS.

# -----------------------------------------------------------------------------
# AMINO ACID METABOLISM
# -----------------------------------------------------------------------------
# O'Neill 2016 highlights amino acids broadly, with specific emphasis on:
# glutamine (anaplerosis, T cell activation), arginine (M1 iNOS vs M2 arginase,
# MDSC immunosuppression), tryptophan (IDO axis, immune tolerance),
# serine/glycine (one-carbon metabolism). This pathway is broad by design.
#
# Biocrates Q1000: Amino acids class (20 proteinogenic), Amino acid-related
# class (including modified / non-proteinogenic AAs: hydroxy-, N-acetyl-,
# methyl-amino acids), Biogenic amines (dopamine, histamine, serotonin, etc. —
# AA-derived), Indoles and derivatives (tryptophan-pathway products),
# Pyridinecarboxylic acids (kynurenine pathway end products: picolinic,
# quinolinic, kynurenic acids), Polyamines (ornithine-derived).
AA_METABOLITE_CLASSES = {
    'Amino acids',
    'Amino acid-related',
    'Biogenic amines',
    'Indoles and derivatives',
    'Pyridinecarboxylic acids',
    'Polyamines',
}
AA_INDIVIDUAL = set()

# =============================================================================
# Final pathway definitions used by the build script
# =============================================================================

PATHWAY_METABOLITES = {
    'IMMUNOMET_ONEIL_GLYCOLYSIS': GLYCOLYSIS_METABOLITES,
    'IMMUNOMET_ONEIL_TCA': TCA_METABOLITES,
    'IMMUNOMET_ONEIL_FAO': set(),  # class-based, see PATHWAY_CLASSES
    'IMMUNOMET_ONEIL_FAS': set(),  # class-based, see PATHWAY_CLASSES
    'IMMUNOMET_ONEIL_AA_METABOLISM': set(),  # class-based
}

PATHWAY_CLASSES = {
    'IMMUNOMET_ONEIL_GLYCOLYSIS': set(),
    'IMMUNOMET_ONEIL_TCA': set(),  # metabolite-based
    'IMMUNOMET_ONEIL_FAO': FAO_METABOLITE_CLASSES,
    'IMMUNOMET_ONEIL_FAS': FAS_METABOLITE_CLASSES,
    'IMMUNOMET_ONEIL_AA_METABOLISM': AA_METABOLITE_CLASSES,
}

# Explicit exclusion of specific metabolites from a pathway even if their class
# would otherwise match. Used to remove biologically incorrect memberships that
# arise from chemical-class grouping (e.g., TCA intermediates sitting in the
# Dicarboxylic acids chemical class).
PATHWAY_EXCLUDE = {
    'IMMUNOMET_ONEIL_GLYCOLYSIS': set(),
    'IMMUNOMET_ONEIL_TCA': set(),
    'IMMUNOMET_ONEIL_FAO': FAO_EXCLUDE,
    'IMMUNOMET_ONEIL_FAS': set(),
    'IMMUNOMET_ONEIL_AA_METABOLISM': set(),
}

PATHWAY_DESCRIPTIONS = {
    'IMMUNOMET_ONEIL_GLYCOLYSIS':
        'Glycolysis: glucose, pyruvate, lactate. Core pathway of activated effector immune cells (Warburg-like). '
        'Framework: O\'Neill, Kishton, Rathmell 2016 Nat Rev Immunol.',
    'IMMUNOMET_ONEIL_TCA':
        'Tricarboxylic acid cycle intermediates including itaconate (ACOD1 product) and 2-hydroxyglutarate. '
        '"Broken" TCA with citrate/succinate accumulation is a hallmark of M1 macrophages. '
        'Framework: O\'Neill et al. 2016.',
    'IMMUNOMET_ONEIL_FAO':
        'Fatty acid oxidation: acylcarnitines (carnitine shuttle products) and dicarboxylic acids '
        '(omega-oxidation products). FAO dominates in memory T cells, M2 macrophages, and Tregs. '
        'Framework: O\'Neill et al. 2016.',
    'IMMUNOMET_ONEIL_FAS':
        'De novo fatty acid synthesis products and related glycerolipids/choline-containing phospholipids. '
        'FAS supports membrane expansion in activated DCs and Th17 differentiation. '
        'Framework: O\'Neill et al. 2016.',
    'IMMUNOMET_ONEIL_AA_METABOLISM':
        'Amino acid metabolism: 20 proteinogenic amino acids, amino acid derivatives, biogenic amines, '
        'kynurenine pathway products, and polyamines. Includes glutamine (T cell activation), arginine '
        '(M1/M2 polarization), tryptophan (IDO-kynurenine axis). '
        'Framework: O\'Neill et al. 2016.',
}

# Order for stable output
PATHWAY_ORDER = [
    'IMMUNOMET_ONEIL_GLYCOLYSIS',
    'IMMUNOMET_ONEIL_TCA',
    'IMMUNOMET_ONEIL_FAO',
    'IMMUNOMET_ONEIL_FAS',
    'IMMUNOMET_ONEIL_AA_METABOLISM',
]
