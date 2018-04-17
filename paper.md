---
title: 'GMSimViz: Automated 3D Visualization of Ground Motion Simulation with Generic Mapping Tools (GMT)'
tags:
  - ground motion simulation
  - 3d visualization
  - automated workflow
authors:
- name: Viktor Polak
   orcid: 0000-0000-0000-0000
   affiliation: "1"
- name: Brendon Bradley
   orcid: 0000-0000-0000-0000
   affiliation: 2
affiliations:
- name: QuakeCoRE, University of Canterbury, Christchurch, New Zealnd
   index: 1
- name: Department of Civil and Natural Resources Engineering, University of Canterbury, Christchurch, New Zealand
   index: 2
date: 22 March 2018
bibliography: paper.bib
---

# Summary
GMSimViz is an automation tool that produces an animated 3D visualization of geological faults, ground motion and other earthquake related data. Typically ground motion simulations are computed by a High Performance Computing (HPC) facility, and its verification involves various data visualization methods.
A 3D animation of the ground motion is an excellent media to understand the nature of an earthquake, and to communicate with the general public. However, its production has been largely left to time-consuming manual interaction with a 3D visualization software package, such as Paraview as no existing solution provides a fully-automated workflow.

GMSimViz was created to provide a fully automated workflow to produce a quality 3D animation directly from the ground motion simulation data. It uses the Generic Mapping Tools (GMT) to create individual frames that are joined by FFmpeg to create a movie file.

GMT is a collection of command line interface (CLI) programs. Creating a fully featured 3D animation with a high frame rate often involves execution of millions of GMT commands. GMSimViz has a Python wrapper to help automate the rendering process through scripting. Various visual effects such as tilt and rotation of the view and fading in/out of different elements are fully automated. The view window and related variables such as viewing angle, center, level of zoom etc. are determined algorithmically based on the data extracted from the input files.

GMSimViz extends the 3D support provided by GMT by using GMT datum to projection coordinate conversion. This allows plotting vertical or near-vertical surfaces in any orientation.
The output animation primarily contains the view of geological faults and ground motion. Users can optionally add map data for liquefaction and landslide probability as well as road network status data.

The rendering process is computationally intensive, yet provides an excellent opportunity to compute in parallel. GMSimViz uses Message Passing Interface (MPI) to utilize a multi-core computer or High Performance Computing (HPC) facilities.
