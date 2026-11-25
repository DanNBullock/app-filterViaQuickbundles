[![Abcdspec-compliant](https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg)](https://github.com/brain-life/abcd-spec)
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-bl.app.###-blue.svg)](https://doi.org/10.25663/bl.app.###)

# [app name]

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

> **Removes spurious streamlines from a .tck tractogram by clustering it with DIPY QuickBundlesX and culling streamlines in small clusters, also splitting an optional white matter classification.**

Clusters the whole-brain tractogram with dipy.segment.bundles.qbx_and_merge and keeps only streamlines in clusters with more than streamThresh members; the rest are written to a separate culled tractogram. If a white matter classification is supplied, its labels are split to match the surviving and culled streamlines. Use it as a quick outlier-removal step after tractography. The README is an unfilled template, and qbThreshes is parsed but not passed to the clustering, which always uses thresholds 30, 20, 10, 5 mm.

| | |
|---|---|
| 📥 **Inputs** | `neuro/track/tck` · `neuro/wmc` (optional) |
| 📤 **Outputs** | `neuro/track/tck` ×2 · `neuro/wmc` ×2 |
| 🐳 **Container** | `brainlife/dipy:1.4.1` |
| ⚙️ **Resources** | 4 CPUs · 30 GB memory · 4 h walltime |
| 🔗 **Run it** | [doi.org/10.25663/brainlife.app.564](https://doi.org/10.25663/brainlife.app.564) |

#### 📥 Inputs

| Input | Datatype | Required | Description |
|---|---|---|---|
| `tractogram` | `neuro/track/tck` | yes | A whole brain tractogram that is to be subjected to a cull of extraneous or otherwise spurious streamlines |
| `classification` | `neuro/wmc` | no | Optional white matter classification (classification.mat) indexing the tractogram's streamlines; its labels are carried over to the surviving and culled outputs. |

#### 📤 Outputs

| Output | Datatype | Description |
|---|---|---|
| `tck_survived` | `neuro/track/tck` | The output tractogram which has had some number of spurious or extraneous streamlines removed. |
| `wmc_survived` | `neuro/wmc` | A WMC corresponding to the surviving streamlines. Identities passed through if WMC input provided, otherwise ignore. |
| `tck_culled` | `neuro/track/tck` | A tractogram containing the streamlines that were culled |
| `wmc_culled` | `neuro/wmc` | A WMC corresponding to the culled streamlines. Identities passed through if WMC input provided, otherwise ignore. |

<details>
<summary><b>⚙️ Parameters</b> (2)</summary>

| Parameter | Type | Default | Description |
|---|---|---|---|
| `streamThresh` | number | `5` | The minimum streamlines (and thus an **integer**) needed in a quickBundles cluster in order to survive the culling process. |
| `qbThreshes` | string | `32,16,8,3` | Comma-separated QuickBundlesX distance thresholds; parsed by the code but not passed to qbx_and_merge, so the defaults 30,20,10,5 are always used. |

</details>

<!-- brainlife-app:end -->

## Overview description

This is where the description of the app goes.  Provide quick summaries of:

### App use case

### Overly simplified algorithm/metholodology

### Necessary inputs and outputs

_Check out the [brainlife datatypes webpage](https://brainlife.io/datatypes) for a cataloguing of relevant datatypes._

## Author, funding sources, references, & license info

### Authors
- [Dan Bullock](https://github.com/DanNBullock/) ([bullo092@umn.edu](mailto:bullo092@umn.edu))

### PI
- [Sarah Heilbronner](https://med.umn.edu/bio/department-of-neuroscience/sarah-heilbronner) ([heilb028@umn.edu](mailto:heilb028@umn.edu))

### Funding
[![NIH-NIBIB-1T32EB031512-01](https://img.shields.io/badge/NIH_NIBIB-1T32EB031512--01-blue.svg)](https://reporter.nih.gov/project-details/10205698)
[![NIMH-5P50MH119569-02](https://img.shields.io/badge/NIMH-5P50MH119569--02-blue.svg)](https://reporter.nih.gov/project-details/10123009)
[![NIMH-5R01MH118257-04](https://img.shields.io/badge/NIMH-5R01MH118257--04-blue.svg)](https://reporter.nih.gov/project-details/10122991)
[![NIDA-1P30DA048742-01A1](https://img.shields.io/badge/NIDA-1P30DA048742--01A1-blue.svg)](https://reporter.nih.gov/project-details/10025457)

* [Dan Bullock](https://github.com/DanNBullock/)'s work is supported by the following sources:
    - The University of Minnesota’s Neuroimaging Postdoctoral Fellowship Program, College of Science and Engineering's, and the Medical School's NIH funded T32 Neuroimaging Training Grant. NOGA: [1T32EB031512-01](https://reporter.nih.gov/project-details/10205698) 

- [Sarah Heilbronner](https://med.umn.edu/bio/department-of-neuroscience/sarah-heilbronner)'s work is supported by the following sources:
    - The University of Minnesota’s Neurosciences' and Medical School's NIMH grant for the investigation of "Neural Basis of Psychopathology, Addictions and Sleep Disorders Study Section[NPAS]". NOGA: [5P50MH119569-02-04](https://reporter.nih.gov/project-details/10123009) 
    - The University of Minnesota’s Neurosciences' Translational Neurophysiology grant. NOGA: [5R01MH118257-04](https://reporter.nih.gov/project-details/10122991)
    - The University of Minnesota’s Neurosciences' Addiction Connectome grant. NOGA: [1P30DA048742-01A1](https://reporter.nih.gov/project-details/10025457)
    
### Citations

#### Directly relevant citations

(Provide citations that are _directly_ relevant to this code implementation here)

#### Indirectly relevant citations

(Provide citations that are _indirectly_ relevant to this code implementation here)

### License info

(Statement of copyright and license)

## Nuts and bolts -- Using this app
(Provide a description of how to use this code repository in different contexts)

### Omnibus usage notes
(provide notes/descriptions/documentation of characteristics and info that is applicable _irrespective_ of the particulars of the compute environment/use-case)

#### The config.json 

The [config.json file](https://brainlife.io/docs/apps/helloworld/#configjson) is a [standard component of application functionality on the brainlife.io platform](https://brainlife.io/docs/apps/introduction/).  However, even outside of the brainlife.io context (running this code in your local python environment), the config.json file is a clean and effective way of managing file and parameter inputs for this code.  Typically, the "main" python file/script (main.py, is the default case) reads the .json file (using the [json](https://docs.python.org/3/library/json.html) module) into a _config_ [dictionary object](https://docs.python.org/3/tutorial/datastructures.html#dictionaries) thusly:

```
with open('config.json') as config_json:
	config = json.load(config_json)
```

Variables and parameters can then be read in from the dictionary using the relevant _keys_.  An example config.json setup might look something like this:

```
{
    "inputFile_1": "local/path/to/inputFile_1",
    "inputFile_2": "local/path/to/inputFile_2",
    "inputFile_3": "local/path/to/inputFile_3",
    "parameter_1": parameter_1_value,
    "parameter_2": parameter_2_value,
    "parameter_3": parameter_3_value
}
```
Consider reviewing the [json standard overview](https://www.json.org/json-en.html) for help formatting this object.

The config.json file can provide a standard interface for controlling execution of the code, whether using [brainlife.io](https://brainlife.io/), [docker](https://docs.docker.com/)/[singularity](https://sylabs.io/guides/2.6/user-guide/index.html), or a local python environment.

###  Using this app on [Brainlife.io](https://brainlife.io/)

#### Input datatypes

(describe and link the input datatypes)

#### Parameter settings

(parameter settings, including default/typical inputs/values, should be documented and provided on the [brainlife App user interface page](https://doi.org/10.25663/bl.app.###).  However, use this section to also duplicate   )

- parameter_1: parameter 1 description
- parameter_2: parameter 1 description
- parameter_3: parameter 1 description

#### Single execution notes

(notes pertinent/specific to usage via the **Execute** interface)

#### Rule/Pipeline usage notes

(notes pertinent/specific to usage via the **Pipeline/Rule** interface)

### LOCAL usage via [docker](https://docs.docker.com/)/[singularity](https://sylabs.io/guides/2.6/user-guide/index.html)

Although excution of python code/apps (under the current developmental framework) is typically controlled by the _main.py_ file, for the purposes of portability and standardization, the execution of of this file is acheived via a bash script that sets up some environmental variables for local / HPC usage and then runs the python script in the relevant docker environment (using "[singularity exec](https://sylabs.io/guides/3.1/user-guide/cli/singularity_exec.html)") in (roughly) the following fashion:

```
#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=02:00:00

# run the actual python code
singularity exec docker://organization/container:version python3 main.py
```

#### Installation
Local use of [docker images/containers](https://www.techtarget.com/searchitoperations/definition/Docker-image) using [singularity](https://sylabs.io/guides/2.6/user-guide/index.html) requires local installation of singularity (a non-trivial matter).  The [Sylabs Singularity documentation page](https://sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps) provides an overview of installation processes, however, this may not cover all installation cases.  For example [a debian package](https://neuro.debian.net/install_pkg.html?p=singularity-container) has been provided by the [neurodebian group](https://neuro.debian.net/).  Additionally, [Brainlife documention](https://brainlife.io/docs/) provides a [guide to singularity installation/usage as well](https://brainlife.io/docs/apps/container/).

#### Settings & setup

There are two major singularity properties to keep in mind when considering local Singularity configuration: the [cache directory](https://sylabs.io/guides/3.3/user-guide/build_env.html) and the [bind path](https://sylabs.io/guides/3.3/user-guide/bind_paths_and_mounts.html)

##### Singularity cache directory
Docker images that are pulled from repositories like [dockerhub](https://hub.docker.com/) are stored and "built" locally.  Usage of a wide range of docker images can quite easily lead to a significant buildup of files in this directory, which can quickly occupy a great deal of harddrive.  As such it is important to be mindful of where (e.g. which disk resource) this environmental variable points to.

##### Singularity bind path
The virtual environments instantiated when a docker image is run (via singularity) **do not** have full access to your local file system.  Instead, only [a narrow set of paths](https://sylabs.io/guides/2.6/user-guide/bind_paths_and_mounts.html#system-defined-bind-points) are made available by default.  As such, **any functions or files called by the code that ARE NOT on these paths will not be found, resulting in an error**.  There are two approaches to dealing with this:

- 1.  Ensure that all requisite files & functions are accessible on these paths.  For requisite files (e.g. input data) in particular, this may mean storing these in a subdirectory of the directory from which the "main" bash script (e.g. where the singularity call is executed).
- 2.  [Manual specification of additional bindpaths using the --bind option](https://sylabs.io/guides/3.1/user-guide/bind_paths_and_mounts.html#user-defined-bind-paths)

### Local usage via python
With the appropriate modules installed in your local python environment it is also possible to run this code.  Ensure that the required modules (see below) are installed, that the config.json file is pointing to the correct files, and you should be ready to go.

#### Required packages / modules
(list the packages and modules required for the execution of this code)

(consider https://docs.python.org/3/library/modulefinder.html run on main)
(non-exhaustive listing of directly called packages that the relevant docker image would need to have or the user would need to have in their local environment)

## Followup 
(consider describing and linking code and applications which could be used with the output of this code or application to produce a more comprehensive analysis or processing pipeline)

## Development note
Work to better structure and document this code/application.  Please feel free to create issues [using the "issues" tab above] to help foster clarity in documentation, or to suggest alterations to documentation/code more directly. Furthermore, it is _acutely_ understood that many of the functionalities in this package may be redundant implementations of (likely more robust and reliable) functionalities from other packages. The identification of such instances (and their other-package correspondances) would be **greatly appreciated**. Feel free to create branches which implement these alternatives.  Be sure to update the documentation as well to describe your changes and to ensure that your contributions are appropriately credited.

