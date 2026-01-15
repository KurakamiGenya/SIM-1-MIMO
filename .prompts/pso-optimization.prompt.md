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

**System to Optimize**: {{SYSTEM_NAME}}
- Problem dimension: {{DIMENSION}} (example: 1000-2000 variables)
- Variable type: {{VARIABLE_TYPE}} (example: phase angles θ ∈ [-π, π])
- Objective: {{OBJECTIVE}} (example: Maximize spectral efficiency)

**Baseline Algorithm**: {{BASELINE_NAME}}
- Current implementation: {{BASELINE_FILE}}
- Performance metrics: {{BASELINE_METRICS}}

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
- **Main implementation**: `{{PREFIX}}_PSO.m` - Full optimization with all layers/scenarios
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

---

## Template Variables

Fill in these before using:
- `{{SYSTEM_NAME}}`: Name of system being optimized
- `{{DIMENSION}}`: Number of variables/dimension
- `{{VARIABLE_TYPE}}`: Type of variables (continuous, discrete, angles, etc.)
- `{{OBJECTIVE}}`: What to maximize or minimize
- `{{BASELINE_NAME}}`: Name of comparison algorithm
- `{{BASELINE_FILE}}`: Filename of baseline implementation
- `{{BASELINE_METRICS}}`: Key performance metrics
- `{{PREFIX}}`: Prefix for generated files
