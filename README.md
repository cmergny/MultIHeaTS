<div align="center">

# MultIHeaTS

MultIHeaTS is a Multi-layered Implicit Heat Transfer Solver. 

It is an implicit numerical model that simulates and predicts the surface temperature in 1D multi-layered planetary surfaces exposed to solar radiation.

[Getting started](#getting-started) •
[Installation](#installation) •
[How to Use](#how-to-use) •
[License](#license)

</div>



# Getting Started

 ![Image description which will be the alt tag](examples/temp_evo.gif)


## Dependancies

- a shell
- git
- conda

You can find conda at https://www.anaconda.com/ although I would suggest installing it directly from the command line.

## Installation

Copy the project localy using git clone:

```bash
git clone git@gitlab.dsi.universite-paris-saclay.fr:cyril.mergny/multiheats.git
```
then cd to the path of the repositery on you computer

```bash
cd path_to_multiheats/
```

Then install the required conda environment:

```bash
conda env create -f environment.yml
```

# Contributing

Contributions are welcome:

- Feel free to open an issue for feedback about usability.
- You may fork the project as you wish as long as you cite the original in your research.
- Pull request may be accepted if new features are in the scope of the MultIHeaTS core.

Please keep pull requests focused and don't change multiple things at the same
time.



# License

MultIHeaTS is distributed under the terms of the GNU GPL License Version 3. A complete version of the license is available in the COPYING file in this repository. Any contribution made to this project will be licensed under the GNU GPL License Version 3.
