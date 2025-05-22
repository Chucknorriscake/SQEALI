# Welcome to SQEALI

Simulation of QUEstionable ALgorythmic Integrity

## Description

This framework streamlines the creation of MEEP simulations, transitioning from a code-based approach to an easily editable configuration file approach, where all simulation parameters are handled through .json files.

## Installation

Prerequisites for this framework include a running version of Anaconda or Miniconda installed on your computer.
To install the environment, run:

```bash
conda create --name <env> --file requirements.txt
conda activate <env>
```

## Prerequisits

The framework was developed and tested in Python 3.11.11. The main focal packages are:

- pymeep (multithreaded)
- numpy
- scipy
- matplotlib

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.

## Project Status

This is a very crude work-in-progress repository. The codebase is, at best, spaghetti and on fire. The proposed workflow has not been thoroughly tested, and edge cases and error handling are not fully implemented. However, as this project is a fun proof of concept and an ode to OOP—as well as sinking more time into a wonderful programming task than necessary—I am perfectly okay with this and happy about it.

## Author 
Julian Verhey