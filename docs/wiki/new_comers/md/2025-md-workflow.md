---
authors: Zhenghao Wu
title: MD Workflow
comments: true
---

# General Workflow for Molecular Dynamics

This tutorial contains hands-on exmaples for:

1. Simple Lennard-Jones liquid 
2. Coarse-grained bead-spring polymer model
3. All-atom model of polymer in water

It provides a step-by-step guide for **running classical molecular dynamics (MD) simulations using LAMMPS**. It is aimed at beginners and covers environment setup, required software installation, molecular modeling, simulation script preparation, job submission, and data analysis. Please follow each step closely and refer to documentation or community resources when needed.

---

## Table of Contents

- [General Workflow for Molecular Dynamics](#general-workflow-for-molecular-dynamics)
  - [Table of Contents](#table-of-contents)
  - [Environment Setup and Software Installation](#environment-setup-and-software-installation)
    - [Set Up Conda Environment](#set-up-conda-environment)
    - [Other Essential Software](#other-essential-software)
  - [Detailed LAMMPS Examples](#detailed-lammps-examples)
    - [Simple Lennard-Jones liquid](#simple-lennard-jones-liquid)
    - [Coarse-grained bead-spring polymer model](#coarse-grained-bead-spring-polymer-model)
    - [All-atom model of polymer in water](#all-atom-model-of-polymer-in-water)
  - [FAQ and Useful Resources](#faq-and-useful-resources)
  - [Notes](#notes)

---

## Environment Setup and Software Installation

### Set Up Conda Environment

- **Module Anaconda/Miniconda:**
    ```bash
    module load miniconda/3
    ```
  - see detail from the [official site]([https://repo.anaconda.com/](https://wugroupwiki.github.io/wiki/cluster_usage/conda/#_1)) and follow our group wiki.
- **Create and activate a Python environment:**

    ```bash
    conda create --name py311 python=3.11
    conda activate py311
    ```

- **Install required packages:**

    ```bash

    # If network issues occur, change sources in ~/.condarc to Tsinghua mirrors, e.g.:
    # channels:
    # - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
    # - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
    # - conda-forge
    # - defaults
    # ssl_verify: true
    ```

### Other Essential Software

- **Download and install:** [`packmol`](http://m3g.iqm.unicamp.br/packmol), [`pubplot`](https://github.com/Chenghao-Wu/pubplot)
- **Tip:** Always read each software's README before installing.

---

## Detailed LAMMPS Examples

This set of tutorials consists of three tutorials arranged in order of increasing difficulty. Although each tutorial can be read independently, information introduced in earlier tutorials is generally not repeated in detail in later ones. For this reason, we recommend that beginners follow the tutorials in order. The novelties associated with each tutorial are briefly described below.

### [Simple Lennard-Jones liquid](lennardjones-md.md)

> **Introduction**: In Lennard-Jones fluid, the structure of LAMMPS input files is illustrated through the creation of a simple atomic Lennard-Jones fluid system. Basic LAMMPS commands are used to set up interactions between atoms, perform an energy minimization, and finally run a simple MD simulation in the microcanonical (NVE) and canonical (NVT) ensembles.

### [Coarse-grained bead-spring polymer model](polymer-md.md)

> **Introduction**: The objective of this tutorial is to perform simple MD simulations using LAMMPS. The system consists of a Lennard-Jones fluid composed of neutral particles with two different effective diameters, contained within a cubic box with periodic boundary conditions. In this tutorial, basic MD simulations in the microcanonical (NVE) and canonical (NVT) ensembles are performed, and basic quantities are calculated, including the potential and kinetic energies.

### [All-atom model of polymer in water](polymer-in-water-md.md)

> **Introduction**: In Polymer in water, two component, liquid water (flexible three-point model) and a polymer molecule, are merged and equilibrated. A long-range solver is used to handle the electrostatic interactions accurately, and the system is equilibrated in the isothermal-isobaric (NPT) ensemble; then, a stretching force is applied to the polymer. Through this relatively complex solvated polymer system, the tutorial demonstrates how to use type labels to make molecule files more generic and easier to manage [3].

## FAQ and Useful Resources

- **Encountering Errors:**  
  Use Google, Computational Chemistry forums, and official documentation for troubleshooting.  
  For LAMMPS-specific errors, consult the [LAMMPS Error Guide](https://docs.lammps.org/Manual.html#errors).

- **More Resources:**
    - [LAMMPS Manual](https://docs.lammps.org/Manual.html)
    - [XJTLU HPC Usage Guide](https://hpc.xjtlu.edu.cn/)
    - [OPLS-AA and SMILES for LAMMPS Tutorial](https://longkunxuluke.github.io/posts/2020/11/blog-post-4/)
    - [Draw structure and get SMILES (Marvin JS)](https://marvinjs-demo.chemaxon.com/latest/demo.html)
    - [VMD - Visual Molecular Dynamics](https://www.ks.uiuc.edu/Research/vmd/)
    - [Bilibili Video with MD theory and practices](https://www.bilibili.com/video/BV1ZWiiYfESi/?spm_id_from=333.1387.homepage.video_card.click&vd_source=33fe2c9e68d4739a445f949401b346e1)

---

## Notes

- Make good use of bash aliases and scripts to accelerate batch operations.
- Always maintain organized file/script naming and archiving habits for reproducibility and collaboration.
- If in doubt at any step, consult the relevant README or Wiki, or ask more experienced colleagues.

---

*This tutorial is subject to updates. Contributions and suggestions are welcome!*
