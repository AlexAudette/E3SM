# Running E3SM with source-aware numerical water tracers

## General information

iE3SM is a version of the Department of Energy Energy Exascale Earth System Model version 2.1 (E3SM v2.1) that provides a source-aware water tagging capability. By defining source regions using the `setup_sim.py` tool, the user can define as many regions as desired for any latitude and longitude bounds. See Audette et al. (2026) for more details. 

[Audette, A., Feldl, N., Singh, H., Heyblom, K., Wang, H., Nusbaumer, J., et al. (2026). Numerical water tracers in the atmospheric component of the Energy Exascale Earth System Model: Implementation and changes in moisture origin. Journal of Advances in Modeling Earth Systems, 18, e2025MS005287.](https://doi.org/10.1029/2025MS005287)

## Getting the code

Clone the water tracer fork of E3SM:

```bash
git clone -b ie3sm1.0 --recursive https://github.com/AlexAudette/E3SM.git ~/E3SM_wtrc
cd E3SM_wtrc
```

If you are on Perlmutter, no further machine configuration is needed. If you are on another machine, it is recommended to first install and validate the `maint-2.1` branch of the standard E3SM to confirm your machine is set up correctly before switching to this fork. See the [General E3SM information](#general-e3sm-information) section for instructions on how to do this.

---

## Creating a new case

### 1. Define your case directory

Set a variable pointing to where your case will live. On Perlmutter, my cases are stored under `/global/cfs/cdirs/<project>/<username>/E3SM_simulations`:

```bash
export CASEDIR="/global/cfs/cdirs/<project>/<username>/E3SM_simulations/<casename>"
```

### 2. Create the case

Run `create_newcase` from the E3SM `cime/scripts` directory:

```bash
/global/homes/<initial>/<username>/E3SM_wtrc/cime/scripts/create_newcase --case $CASEDIR --compset 1850SOI_EAM%CMIP6_ELM%SPBC_MPASSI_MPASO_MOSART_SGLC_SWAV --res ne30pg2_EC30to60E2r2 --mach pm-cpu --project m4426
```

The key options are:

| Option | Value | Description |
|---|---|---|
| `--case` | `$CASEDIR` | Path where the case directory will be created |
| `--compset` | `1850SOI_EAM%CMIP6_...` | Pre-industrial control compset with the EAM atmosphere and MPAS ocean/ice |
| `--res` | `ne30pg2_EC30to60E2r2` | Atmosphere at ~1° (ne30pg2) with the EC30to60 MPAS ocean mesh |
| `--mach` | `pm-cpu` | Perlmutter CPU partition |
| `--project` | `m4426` | Compute allocation to charge |

### 3. Configure the water tracers with `setup_sim.py`

The `setup_sim.py` script is a Python tool that automates all water tracer configuration. It can be cloned from my github repository:

```bash
git clone https://github.com/AlexAudette/iE3SM_utils.git
```

It should be ran form within the `iE3SM_utils/` directory — it writes directly into `$CASEDIR`, which is passed as a command-line argument.

The `iE3SM_utils/` directory contains the Python tooling that automates all water tracer configuration. It should be ran form within the `iE3SM_utils/` directory — it writes directly into `$CASEDIR`, which is passed as a command-line argument.

#### What it does

Running `setup_sim.py` performs three steps automatically:

1. **Writes `user_nl_eam`** — generates the EAM namelist with tracer names, isotope ratios, and output field lists (`fincl1`/`fincl2`), writing it to `$CASEDIR/user_nl_eam`.
2. **Patches `atm_import_export.F90`** — inserts the tracer source definition block (geographic bounds and flux assignments) into `$CASEDIR/SourceMods/src.eam/atm_import_export.F90`.
3. **Runs `xmlchange` commands** — configures the run length, resubmission count, wall-clock time, and archiving settings for the case.

#### Directory structure


```
iE3SM_utils/                          ← run scripts from here
├── setup_sim.py
├── tracer_config.py
├── namelist.py
├── f90_patch.py
├── xml_config.py
├── tracer_configuration.json       ← edit this to define your tracers
└── run_config.json                 ← edit this to set run length, wallclock, etc.

$CASEDIR/                           ← files will be written here
├── SourceMods/
│   └── src.eam/
│       └── atm_import_export.F90   ← will be patched by setup_sim.py
└── user_nl_eam                     ← will be written by setup_sim.py
```

#### Configuration

Before running, review and edit the two JSON files:

**`tracer_configuration.json`** — defines the water tracer regions. `H2O` must always be the first entry. Add or modify regional tags (e.g. `ARC`, `ROW`) with their latitude/longitude bounds and isotope ratios. Here are examples for the two tracers describing the moisture evaporated from the Arctic and from the rest of the world:

```json
{
    "H2O": {
        "lat_bounds": [-90, 90],
        "rstd": 1.0,
        "description": "Bulk water reference. Required in every water tracer run."
    },
    "ARC": {
        "lat_bounds": [70, 90],
        "rstd": 0.003,
        "description": "Arctic water vapour"
    },
    "ROW": {
        "lat_bounds": [-90, 70],
        "rstd": 0.997,
        "description": "Rest of world water vapour"
    }
}
```

**`run_config.json`** — controls the CESM case settings:

```json
{
    "resubmit": 5,
    "stop_n": 5,
    "stop_option": "nyears",
    "wallclock": "48:00:00",
    "dout_s": true
}
```

#### Configure the water tracer files

From the `iE3SM_utils/` directory, pass `$CASEDIR` as an argument:

```bash
python setup_sim.py $CASEDIR
```

The script will print a confirmation line for each step as it completes:

```
Written user_nl_eam
Patched SourceMods/src.eam/atm_import_export.F90 successfully.
Running: ./xmlchange --append CAM_CONFIG_OPTS=-water_tracer h2o -water_tag_num <N_WTRC>
Running: ./xmlchange RESUBMIT=5
...
```

Once this has completed without errors, the case is fully configured and ready to set up and build:

### 4. Set up and build the case

```bash
cd "${CASEDIR}"
```
At this point you might want to adjust the `$CASEDIR/env_mach_pes.xml` which contains the PE layout for the run. You can always change it later, but you will need to recall the `case.setup` script with the `--reset` option, which will reset the `JOB_WALLCLOCK_TIME` variable in `env_workflow.xml`. 


```bash
./case.setup

./case.build
```

At this point, the model should have built sucessfully and the console should show a message indicating so. The case is now ready to be submitted to the queue or ran interactively on an interactive compute node (not recommended). 


### 5. Submit the case

```bash
./case.submit
```

This will submit the job to the queue.

---

## General E3SM information

### Getting the standard E3SM code

```bash
git clone -b maint-2.1 --recursive https://github.com/E3SM-Project/E3SM.git
cd E3SM
```

Available branches include `maint-1.0`, `maint-1.1`, `maint-1.2`, `maint-2.0`, and `maint-2.1`. Do not use `.tar.gz` or `.zip` archives from GitHub — they do not include the submodule code. `maint-2.1` is the release from which this project is branched from. Submodules such as CIME have not been modified in iE3SM, and are maintained by the E3SM team. 

### Input data

On supported platforms (including Perlmutter), all required input data is already available in world-readable local directories. On unsupported machines, do not attempt to download all input data at once — the CIME case control system will automatically fetch any missing files from the E3SM input data server (`https://web.lcrc.anl.gov/public/e3sm/`) or the CESM input data server as each case requires them.

### Further reading

- [E3SM Step-by-Step running guide](https://docs.e3sm.org/running-e3sm-guide/)
- [E3SM quick start page](https://e3sm.org/model/running-e3sm/e3sm-quick-start/)

### Troubleshooting

- Report a bug: https://github.com/E3SM-Project/E3SM/issues
- Ask a question: https://github.com/E3SM-Project/E3SM/discussions

---
MIT License

Copyright (c) 2026 Alexandre Audette

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.