---
description: Implement Standard PSO for optimization problems with gradient-based baseline comparison
author: duc
created: 2026-01-12
tags: [optimization, PSO, particle-swarm, gradient-descent, comparison]
---

# Standard PSO Implementation and Comparison

You are an expert in optimization algorithms and MATLAB programming. Help me implement Standard Particle Swarm Optimization (PSO) for a complex optimization problem.

## Problem Context

I have an optimization problem with the following characteristics:

**System to Optimize**: Stacked Intelligent Metasurfaces (SIM) for MIMO precoding
- Problem dimension: 1100-2000 variables (100×L + 100×K, where L=1-10 TX layers, K=10 RX layers)
- Variable type: Phase angles θ ∈ [-π, π] for metasurface elements
- Objective: Maximize spectral efficiency (Sum Rate in bits/s/Hz)

**Baseline Algorithm**: Gradient Descent with Random Initialization
- Current implementation: test_SIM.m
- Performance metrics: NMSE (Normalized Mean Square Error), Capacity (Sum Rate), Convergence speed

## Implementation Requirements

### 1. PSO Mapping
Map my problem to Standard PSO framework:
- **Particle representation**: Define how particles encode solution
- **Fitness function**: Define objective to maximize/minimize
- **Constraints**: Handle variable bounds and constraints

### 2. Standard PSO Formula
Use the canonical PSO update equations:
```
v = w*v + c1*rand*(pBest - x) + c2*rand*(gBest - x)
x = x + v
```

Where:
- w = 0.7 (CONSTANT - Standard PSO characteristic)
- c1 = 1.5 (Cognitive parameter)
- c2 = 1.5 (Social parameter)

### 3. Code Structure
Create these files:
- **Main implementation**: `test_SIM_PSO.m` - Full optimization with all layers/scenarios
- **Quick demo**: `demo_PSO.m` - Faster version for testing
- **Comparison script**: `compare_results.m` - Compare PSO vs baseline
- **Documentation**: `PSO_README.md` - Detailed explanation

### 4. Comparison Analysis
Compare PSO with baseline algorithm on:
- Performance metrics (capacity, accuracy, error, etc.)
- Convergence behavior (smooth vs noisy)
- Computational cost (time, iterations, evaluations)
- Scalability with problem dimension
- When each algorithm performs better

## Expected Deliverables

### Code Files
1. ✅ Main PSO implementation with proper parameters
2. ✅ Demo version for quick testing
3. ✅ Comparison and visualization scripts
4. ✅ Documentation with theory and usage

### Analysis
1. ✅ Performance comparison tables
2. ✅ Plots showing algorithm behavior
3. ✅ Discussion of why one algorithm outperforms
4. ✅ Recommendations for practical use

### Issues to Handle
- Fix any toolbox dependencies (use standard MATLAB functions)
- Ensure fitness function matches baseline calculation
- Handle numerical stability (NaN, Inf checks)
- Proper constraint handling (bounds, wrapping)

## Key Questions to Address

1. **Mapping**: How to encode my specific problem into PSO particles?
2. **Fitness**: What is the correct objective function calculation?
3. **Comparison**: Is PSO better or worse than baseline? Why?
4. **Practical**: Which algorithm should I use in practice?
5. **Improvements**: How can PSO performance be improved if needed?

## Success Criteria

- ✅ PSO correctly implements Standard PSO (w=0.7 constant)
- ✅ Fitness function matches baseline methodology
- ✅ Comprehensive comparison showing clear winner
- ✅ Understanding of why one algorithm wins
- ✅ Reproducible results with clear documentation

## Style Preferences

- Explain concepts clearly in Vietnamese when discussing
- Use English for code and technical terms
- Provide concrete examples and visualizations
- Fix issues immediately when they arise
- Give clear step-by-step instructions for running code
