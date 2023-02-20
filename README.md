---
title: "README"
author: "Eloy"
date: "1/2/2023"
output: html_document
---

## SIMULATING TUBULAR TISSUES IN 3D

The code in this repository is used to perform 3D simulations of cellular tissues with tubular geometry. This involves the simulation of individual cells and the overall tissue structure.

To achieve this, we use Voronoi tessellations to simulate individual cells in the tissue. Voronoi cells divide a space into regions, such that each region contains all points that are closer to one and only one specific cell center.

Additionally, the code uses the Metropolis algorithm to take into account the interactions between cells. This algorithm drives the system towards a thermodynamic equilibrium. It involves randomly proposing changes to the system and accepting or rejecting these changes based on the acceptance criteria (this criteria is based on the cell energy, which depends directly of the area and perimeter of each Voronoi cell). By using Metropolis algorithm, the code is able to achieve an accurate and realistic simulation of cellular tissues.

## Structure of the repository

The repository contains two main algorithms for simulation purposes. The first algorithm is the normal algorithm, which can be found in the "N-cylinder" folder. The second algorithm includes a Bending process and can be found in the "Bending algorithm" folder.

Both algorithms have two versions: a normal version and a parallelized version. The parallelized version allows for simulations to be executed in a distributed manner, which can improve the computational efficiency of the simulations.

In addition to the algorithms, the repository also contains scripts in the "Analysis" folder that are used to analyze the results of the code and produce output that is readable by other programming languages. The "Data" folder contains some sample data that has been executed and can be used for testing purposes.

Overall, this repository is designed to provide a flexible and scalable platform for simulation and analysis of various processes.
