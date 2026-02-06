configfile: "config/config.yml"
NEX_IN = config["nexgddp_in"]
OUT = config["nexgddp_out"]

ERA5_IN = config["era5_in"]
ERA5_OUTDIR = config["era5_outdir"]
ERA5_PACK = config["era5_pack"]

MODELS = config["models"]

SCEN_WINDOWS = {
    "historical": ("1991", "2014"),
    "ssp126": ("2015", "2100"),
    "ssp245": ("2020", "2100"),
    "ssp585": ("2020", "2100"),
}

FUTURE_WINDOWS = config["future_windows"]
SSPS_FUTURE = config["ssps_for_future"]

def scen_file(model, scen):
    y0, y1 = SCEN_WINDOWS[scen]
    return f"{OUT}/{model}/{scen}/HEATSTRESS_mon_{model}_{scen}_{y0}-{y1}.nc"

def future_pack_file(model, scen, win_name):
    return f"{OUT}/{model}/packs/HEATSTRESS_pack_{model}_{scen}_{win_name}.nc"

rule all:
    input:
        ERA5_PACK,
        expand(scen_file("{model}", "historical"), model=MODELS),
        expand(scen_file("{model}", "ssp126"), model=MODELS),
        expand(scen_file("{model}", "ssp245"), model=MODELS),
        expand(scen_file("{model}", "ssp585"), model=MODELS),
        expand(f"{OUT}/{{model}}/packs/HEATSTRESS_pack_{{model}}_pseudoHist_1991-2020.nc", model=MODELS),
        expand(future_pack_file("{model}", "{scen}", "{win}"),
               model=MODELS,
               scen=SSPS_FUTURE,
               win=[w["name"] for w in FUTURE_WINDOWS])

rule era5_monthly_1991_2020:
    output:
        ERA5_PACK
    shell:
        "bash pipeline/heatstress_era5_monthly_snake.sh {ERA5_IN} {ERA5_OUTDIR}"

rule compute_monthly:
    output:
        lambda wc: scen_file(wc.model, wc.scen)
    params:
        model="{model}",
        scen="{scen}"
    shell:
        "bash pipeline/heatstress_nexgddp_monthly_snake.sh {NEX_IN} {OUT} {params.model} {params.scen}"

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
        lambda wc: future_pack_file(wc.model, wc.scen, wc.win)
    params:
        start=lambda wc: next(w["start"] for w in FUTURE_WINDOWS if w["name"] == wc.win),
        end=lambda wc: next(w["end"] for w in FUTURE_WINDOWS if w["name"] == wc.win)
    shell:
        "mkdir -p $(dirname {output}) && "
        "cdo -L -z zip_5 seldate,{params.start},{params.end} {input} {output}"
