# Airborne Space Launcher Conceptual Design

## Overview
This project presents the conceptual design of an airborne launch system for deploying payloads into Low Earth Orbit (LEO), with a strong focus on first-stage reusability.

The system leverages an air-launch architecture to reduce atmospheric losses and improve performance, while integrating recovery strategies to enhance cost efficiency and sustainability.

The work combines multidisciplinary analysis including propulsion, trajectory optimization, aerodynamics, structures, and system-level trade-offs.

## Mission Context
The launcher is designed to:
- Deliver payloads (~250 kg) to LEO (400 km SSO)
- Operate from an airborne platform (e.g., modified cargo aircraft)
- Enable recovery and reuse of the first stage

Key drivers:
- Cost reduction through reusability
- Operational flexibility (air-launch)
- Competitive performance vs ground-based systems

## System Architecture

### Air-Launch Concept
- Carrier aircraft: Boeing 767-class platform
- Release altitude: ~10–11 km
- Reduced drag and gravity losses compared to ground launch
- Launch operations performed over ocean for safety

### Launcher Configuration
- Two-stage-to-orbit (2STO) architecture selected based on ROI analysis
- LOX/RP-1 propulsion system
- Common engine design across stages with different nozzle geometries

Key characteristics:
- First stage: reusable, multi-engine configuration
- Second stage: vacuum-optimized, single engine

### First Stage Reusability
Recovery strategy:
- Boost-back maneuver using main engines
- Controlled re-entry with load constraints
- Parachute-assisted splashdown

Design constraints include:
- Maximum structural loads (< 6 g)
- Controlled velocity for parachute deployment
- Safe touchdown conditions

## Methods

### System-Level Design
- House of Quality (HoQ) used to derive requirements
- Functional decomposition and iterative design process
- Trade-offs between:
  - 2STO vs 3STO
  - Reusability vs cost
  - Performance vs complexity

### Propulsion Design
- Custom LOX/RP-1 rocket engine
- Chamber pressure: 160 bar
- Optimization of:
  - O/F ratio via CEA
  - Nozzle expansion ratio
- Rao nozzle geometry for improved performance

Preliminary regenerative cooling analysis ensures thermal feasibility.

### Trajectory Optimization
- 2-DOF point-mass dynamic model
- Optimization variables:
  - Burn times
  - Flight path angle evolution
  - Initial conditions post-release

Objectives:
- Minimize losses (gravity, drag)
- Satisfy orbital insertion constraints
- Ensure recoverability of first stage

### Aerodynamics and Structures
- Simplified aerodynamic model (CD ≈ 0.2 initial assumption)
- Structural sizing based on load cases:
  - Max dynamic pressure (~70 kPa)
  - Re-entry conditions
- Material selection and wall thickness estimation

### Recovery Modeling
- Boost-back trajectory design
- Parachute sizing and deployment sequence
- Trade-off between propellant allocation and recovery feasibility

### Sensitivity Analysis
- Monte Carlo analysis on uncertainties
- Evaluation of orbit insertion accuracy
- Assessment of robustness to parameter variations

## Results

### Performance
- Successful orbit insertion within acceptable error margins
- Total ΔV requirement ~10.5 km/s (including recovery)

### Reusability Trade-offs
- First-stage recovery significantly impacts mass budget
- Additional propellant required for boost-back (~3.7 tons)
- Trade-off between payload capacity and recoverability

### Cost Analysis
- 2STO configuration provides बेहतर ROI compared to 3STO
- Reusability introduces:
  - Higher initial cost (penalty factor)
  - Lower cost per launch over multiple missions

### Critical Findings
- Recovery is the most constraining subsystem
- Air-launch provides performance benefits but adds integration complexity
- Strong coupling between trajectory design and recovery feasibility

## Implementation
The project includes:
- Modular Matlab-based framework for subsystem analysis
- Integrated system-level iterative design loop
- Trajectory propagation and optimization
- Engine preliminary design and sizing
- Cost and performance trade-off analysis

## Key Concepts
- Air-launch systems
- Reusable launch vehicles
- Trajectory optimization
- ΔV budgeting
- Boost-back recovery
- Rocket propulsion (LOX/RP-1)
- System-level design and trade-offs

## Author
Matteo Portantiolo  
MSc Space Engineering – GNC
