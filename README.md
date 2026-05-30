# Radiotherapy Beam Angle Optimization

![MATLAB](https://img.shields.io/badge/MATLAB-R2023a+-blue.svg)
![Optimization](https://img.shields.io/badge/Optimization-Memetic%20Algorithm-brightgreen.svg)
![Course](https://img.shields.io/badge/Course-Optimization%20Methods-orange.svg)

This repository contains the codebase for a Biomedical Engineering (BME) optimization project developed for the Optimization Methods course at the University of Basel. 

The project solves a highly non-convex clinical problem: finding the optimal intersecting angles for two radiotherapy beams to deliver a lethal dose to a central Tumor while completely sparing an adjacent Organ at Risk (OAR).

## 🧬 The Clinical Problem
In external beam radiotherapy, a single high-intensity radiation beam will destroy all healthy tissue in its path. To solve this, oncologists use multiple intersecting beams of lower intensity. The goal is to optimize the angles ($\theta_1$ and $\theta_2$) so that the intersection (the "crossfire") maximizes the dose over the tumor and minimizes the dose to critical adjacent organs.

## 🧮 Mathematical Formulation
The physical space is abstracted into a 2D 100x100 grid. The radiation beams are modeled analytically using a Gaussian decay function based on the perpendicular distance $d$ to the beam's center line:

$$D(x,y) = I \cdot \exp\left(-\frac{d^2}{2\sigma^2}\right)$$

The Objective Function calculates the Sum of Squared Errors (SSE) across three regional masks, heavily penalizing any radiation that exceeds safe limits in the Organ at Risk ($W_{oar}$):

$$J = J_{tumor} + W_{oar} \cdot J_{oar} + J_{normal}$$

## ⚙️ Optimization Strategy (Memetic Algorithm)
Due to strict OAR constraints, the mathematical landscape is highly non-convex, meaning standard gradient descent algorithms frequently get trapped in local minima. To solve this, the project implements a two-phase hybrid "Memetic" algorithm:

1. **Phase 1: Global Search (Particle Swarm Optimization - PSO)**
   - A swarm of 15 particles explores the bounded $[0, \pi]$ solution space.
   - Escapes local minima and identifies the mathematical "safe valley."
2. **Phase 2: Local Refinement (Nelder-Mead Downhill Simplex)**
   - Initialized exactly around the PSO's best global coordinate.
   - Uses a derivative-free geometric tumbling method (reflection, expansion, contraction) to pinpoint the exact decimal optimal angles.

## 📂 Repository Structure
The project is built entirely in Base MATLAB (no external Optimization or Image Processing toolboxes required). 

- `RadiotherapyCost.m` : The Forward Model and Objective Cost Function.
- `RadiotherapyCostMap.m` : Vectorized wrapper for generating 3D surface visualizations.
- `run_PSO.m` : The bounded Particle Swarm Optimization algorithm.
- `main.m` : The master script running the Memetic pipeline, Nelder-Mead loop, and generating all visualizations.

## 🚀 Getting Started
1. Clone this repository or download the ZIP file.
2. Extract all four `.m` files into a single directory.
3. Open MATLAB and navigate to that directory.
4. Open `main.m` and press **Run**.

## 📊 Visual Outputs
The `main.m` script automatically generates four distinct layers of evaluation:
1. **3D Landscape Tumbling:** Real-time `surf` plotting of the Nelder-Mead simplex navigating the cost function landscape.
2. **Dose Heatmap:** A 2D contour plot showing the physical geometry of the optimized beam crossfire against the organ boundaries.
3. **Clinical DVH:** A standard Dose-Volume Histogram proving 100% target coverage and 0% OAR toxicity.
4. **Pareto Trade-off Front:** A sensitivity analysis looping through various $W_{oar}$ parameter weights to map the fundamental mathematical compromise between tumor underdosing and OAR overdosing.

---
*Developed for the MS in Biomedical Engineering program.*
