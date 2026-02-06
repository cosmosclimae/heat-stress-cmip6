configfile: "config/config.yml"

NEX_IN = config["nexgddp_in"]
OUT = config["nexgddp_out"]

ERA5_IN = config["era5_in"]
ERA5_OUTDIR = config["era5_outdir"]
ERA5_PACK = config["era5_pack"]

MODELS = config["models"]

FUTURE_WINDOWS = config["future_windows"]
SSPS_FUTURE = config["ssps_for_future"]

def scen_file(model, scen):
    if scen == "historical":
        return f"{OUT}/{model}/historical/HEATSTRESS_mon_{model}_historical_1991-2014.nc"
    if scen == "ssp126":
        return f"{OUT}/{model}/ssp126/HEATSTRESS_mon_{model}_ssp126_2015-2100.nc"
    if scen == "ssp245":
        return f"{OUT}/{model}/ssp245/HEATSTRESS_mon_{model}_ssp245_2020-2100.nc"
    if scen == "ssp585":
        return f"{OUT}/{model}/ssp585/HEATSTRESS_mon_{model}_ssp585_2020-2100.nc"
    raise ValueError(f"Unknown scenario: {scen}")

rule all:
    input:
        ERA5_PACK,
        expand(scen_file("{model}", "historical"), model=MODELS),
        expand(scen_file("{model}", "ssp126"), model=MODELS),
        expand(scen_file("{model}", "ssp245"), model=MODELS),
        expand(scen_file("{model}", "ssp585"), model=MODELS),
        expand(f"{OUT}/{{model}}/packs/HEATSTRESS_pack_{{model}}_pseudoHist_1991-2020.nc", model=MODELS),
        expand(f"{OUT}/{{model}}/packs/HEATSTRESS_pack_{{model}}_{{scen}}_{{win}}.nc",
               model=MODELS,
               scen=SSPS_FUTURE,
               win=[w["name"] for w in FUTURE_WINDOWS])

rule era5_monthly_1991_2020:
    output:
        ERA5_PACK
    shell:
        "bash pipeline/heatstress_era5_monthly_snake.sh {ERA5_IN} {ERA5_OUTDIR}"

rule compute_historical:
    output:
        f"{OUT}/{{model}}/historical/HEATSTRESS_mon_{{model}}_historical_1991-2014.nc"
    shell:
        "bash pipeline/heatstress_nexgddp_monthly_snake.sh {NEX_IN} {OUT} {wildcards.model} historical"

rule compute_ssp126:
    output:
        f"{OUT}/{{model}}/ssp126/HEATSTRESS_mon_{{model}}_ssp126_2015-2100.nc"
    shell:
        "bash pipeline/heatstress_nexgddp_monthly_snake.sh {NEX_IN} {OUT} {wildcards.model} ssp126"

rule compute_ssp245:
    output:
        f"{OUT}/{{model}}/ssp245/HEATSTRESS_mon_{{model}}_ssp245_2020-2100.nc"
    shell:
        "bash pipeline/heatstress_nexgddp_monthly_snake.sh {NEX_IN} {OUT} {wildcards.model} ssp245"

rule compute_ssp585:
    output:
        f"{OUT}/{{model}}/ssp585/HEATSTRESS_mon_{{model}}_ssp585_2020-2100.nc"
    shell:
        "bash pipeline/heatstress_nexgddp_monthly_snake.sh {NEX_IN} {OUT} {wildcards.model} ssp585"

rule pseudo_hist_pack_1991_2020:
    input:
        hist=lambda wc: scen_file(wc.model, "historical"),
        ssp126=lambda wc: scen_file(wc.model, "ssp126")
    output:
        f"{OUT}/{{model}}/packs/HEATSTRESS_pack_{{model}}_pseudoHist_1991-2020.nc"
    shell:
        "mkdir -p $(dirname {output}) && "
        "cdo -L -z zip_5 seldate,1991-01-01,2020-12-31 -mergetime {input.hist} {input.ssp126} {output}"

rule future_pack:
    input:
        lambda wc: scen_file(wc.model, wc.scen)
    output:
        f"{OUT}/{{model}}/packs/HEATSTRESS_pack_{{model}}_{{scen}}_{{win}}.nc"
    params:
        start=lambda wc: next(w["start"] for w in FUTURE_WINDOWS if w["name"] == wc.win),
        end=lambda wc: next(w["end"] for w in FUTURE_WINDOWS if w["name"] == wc.win)
    wildcard_constraints:
        scen="ssp126|ssp245|ssp585"
    shell:
        "mkdir -p $(dirname {output}) && "
        "cdo -L -z zip_5 seldate,{params.start},{params.end} {input} {output}"
