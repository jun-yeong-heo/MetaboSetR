# Internal data documentation.
#
# The objects below are stored in `R/sysdata.rda` and built by
# `data-raw/make_sysdata.R`. They are not user-visible. See DESIGN.md
# section 5.2 for the schema.
#
# Objects:
#   metabolite_dict     data.frame  — 1234 metabolite metadata + id annotations
#   status_color_map    data.frame  — hex -> status text mapping
#   pathway_sets        nested list — domain -> set_name -> members
#   pathway_sets_meta   data.frame  — long-format set/member metadata
