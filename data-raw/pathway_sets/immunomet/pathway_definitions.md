# Immunometabolism Pathway Definitions

**Framework:** O'Neill, Kishton, Rathmell (2016). *A guide to immunometabolism
for immunologists*. Nat Rev Immunol 16(9):553–565.

**Supporting references used for individual pathway scoping:**

- Pearce EL, Pearce EJ (2013). Metabolic pathways in immune cell activation
  and quiescence. *Immunity* 38(4):633–643.
- Buck MD, O'Sullivan D, Pearce EL (2015). T cell metabolism drives immunity.
  *J Exp Med* 212(9):1345–1360.
- Mogilenko DA, Sergushichev A, Artyomov MN (2023). Systems immunology
  approaches to metabolism. *Annu Rev Immunol* 41:317–342.

## Scope

Five of the six core pathways in the O'Neill 2016 framework are implemented.
The pentose phosphate pathway (PPP) is excluded because the Biocrates Quant
1000 panel does not cover its canonical metabolites (G6P, 6-phosphogluconate,
ribose-5-phosphate, NADPH).

## Sets

### IMMUNOMET_ONEIL_GLYCOLYSIS

**Canonical metabolites:** Glucose, Pyruvic acid, Lac (lactic acid).

**Coverage on panel:** 3 metabolites.

**Rationale:** O'Neill 2016 Figure 1 defines glycolysis by the substrates and
products at its endpoints. Activated effector immune cells (M1 macrophages,
effector T cells) exhibit aerobic glycolysis with enhanced lactate output.
Panel upstream intermediates (G6P, F1,6BP, PEP) are absent, so the set
captures net flux readouts rather than pathway flux distribution.

### IMMUNOMET_ONEIL_TCA

**Canonical metabolites:** Citric acid, Isocitric acid, AconAcid (aconitate),
a-Ketoglutaric acid, Suc (succinate), Fumaric acid, Malic acid, Oxaloacetic
acid, 2-OH-Glutaric acid, Itaconic acid.

**Coverage on panel:** 10 metabolites.

**Rationale:** TCA intermediates plus the two immunometabolism-relevant
off-cycle metabolites: itaconate (ACOD1/IRG1 product, anti-inflammatory via
SDH inhibition and Nrf2 activation) and 2-hydroxyglutarate (oncometabolite
relevant to IDH-mutant cancers; also produced under hypoxia by wild-type
enzymes). The "broken TCA" state with citrate and succinate accumulation is
a hallmark of M1-polarized macrophages.

### IMMUNOMET_ONEIL_FAO

**Canonical metabolite classes:** Acylcarnitines, Dicarboxylic acids.

**Explicit exclusions from the Dicarboxylic acids class:** 2-OH-Glutaric acid,
a-Ketoglutaric acid, Fumaric acid, Itaconic acid, Malic acid, Oxaloacetic
acid, Suc. These are TCA-specific metabolites that happen to share the
Dicarboxylic acids chemical class but are not products or substrates of
fatty acid β- or ω-oxidation. Biologically they belong to TCA only.

**Coverage on panel:** 58 metabolites (40 acylcarnitines + 18 dicarboxylic
acids after exclusions).

**Rationale:** Acylcarnitines are the diagnostic signature of carnitine-shuttle
flux into mitochondrial FAO. Dicarboxylic acids (excluding TCA intermediates)
include peroxisomal ω-oxidation products (adipic, suberic, DiCA 12:0,
DiCA 14:0) and BCAA-catabolism byproducts that functionally pair with FAO.
FAO dominates in memory T cells, M2 macrophages, and Tregs.

### IMMUNOMET_ONEIL_FAS

**Canonical metabolite classes:** Fatty acids, Triglycerides, Diglycerides,
Monoglycerides, Cholesteryl esters, Phosphatidylcholines,
Lysophosphatidylcholines.

**Coverage on panel:** 445 metabolites.

**Rationale:** O'Neill 2016 treats de novo fatty acid synthesis as the pathway
supporting membrane expansion in activated DCs and Th17 differentiation
(ACC1 dependency). In plasma metabolomics, de novo lipogenesis cannot be
cleanly separated from dietary and lipoprotein-derived pools, so this set is
deliberately broad and captures all plasma lipid readouts that would be
influenced by FAS flux. Users who need a narrower de novo lipogenesis
signature should restrict to SCD-1 substrate/product pairs and related
de-novo-specific metabolites.

**Note on overlap with FAO:** The Fatty acids class is in FAS (as a product).
Acylcarnitines and dicarboxylic acids are in FAO (as FAO-specific readouts).
Individual free fatty acids in plasma reflect net balance of synthesis,
storage, and oxidation; they are assigned to FAS only, which is the
convention adopted here.

### IMMUNOMET_ONEIL_AA_METABOLISM

**Canonical metabolite classes:** Amino acids, Amino acid-related, Biogenic
amines, Indoles and derivatives, Pyridinecarboxylic acids, Polyamines.

**Coverage on panel:** 138 metabolites.

**Rationale:** O'Neill 2016 treats amino acid metabolism as a single broad
pathway with emphasis on glutamine (T cell activation), arginine (M1/M2
macrophage polarization via iNOS vs arginase), tryptophan (IDO-kynurenine
axis for immune tolerance), and serine/glycine (one-carbon metabolism).
The classes included cover all 20 proteinogenic amino acids, modified and
non-proteinogenic derivatives, AA-derived biogenic amines (dopamine,
histamine, serotonin), indole products of the tryptophan pathway, kynurenine
pathway endpoints (picolinic, quinolinic, kynurenic acids), and
ornithine-derived polyamines.

## Overlap handling

Overlap between sets is allowed when biologically meaningful (permissive).
A metabolite may appear in multiple sets when it genuinely participates in
multiple pathways. Artifactual overlap arising from chemical-class grouping
alone (e.g., TCA intermediates appearing in the Dicarboxylic acids class) is
explicitly removed.

**Metabolite-level:** 0 overlaps in the final output — the explicit
TCA-intermediate exclusion from FAO resolves the only artifactual case.

## Direction convention

The immunomet master TSV has no direction column. Immunometabolism does not
have a single reference state (immune cell activation is bidirectional:
effector vs memory, M1 vs M2, Teff vs Treg). Users applying this set in
enrichment analysis should interpret directionality in the context of their
specific comparison
(e.g., activated vs resting T cells, M1 vs M2 macrophages).
