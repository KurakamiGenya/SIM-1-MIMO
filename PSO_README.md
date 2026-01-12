# Standard PSO Implementation for SIM-MIMO

## Mô tả
Thuật toán PSO (Particle Swarm Optimization) chuẩn được áp dụng để tối ưu hóa các góc pha của Stacked Intelligent Metasurfaces (SIM) trong hệ thống MIMO.

## Cấu trúc bài toán

### 1. Input (Hạt - Particle)
Mỗi hạt trong bầy đàn PSO đại diện cho một vector chứa:
- **Góc pha TX-SIM**: M×L góc pha θ cho các phần tử trên các lớp phát
- **Góc pha RX-SIM**: N×K góc pha θ cho các phần tử trên các lớp thu
- **Tổng số chiều**: dim = M×L + N×K = 100×L + 100×10

Mỗi góc pha nằm trong khoảng: θ ∈ [-π, π]

### 2. Hàm mục tiêu (Fitness Function)
**Spectral Efficiency (Sum Rate)** được tính từ công thức:

```
Capacity = Σ log₂(1 + SINR_i)
```

Trong đó:
- SINR_i: Signal-to-Interference-plus-Noise Ratio của stream thứ i
- PA_WF: Power allocation theo thuật toán Water-Filling
- H_SIM: End-to-end channel matrix sau khi qua TX-SIM và RX-SIM

**Mục tiêu**: Maximize Spectral Efficiency (Sum Rate)

### 3. Thuật toán PSO chuẩn

#### Công thức cập nhật vận tốc (Standard PSO):
```matlab
w = 0.7;  % Inertia weight (HẰNG SỐ)
v = w*v + c1*rand*(pBest - x) + c2*rand*(gBest - x);
```

#### Tham số PSO:
- **w = 0.7**: Hệ số quán tính (constant - đặc điểm của PSO gốc)
- **c1 = 1.5**: Hệ số học tập cá nhân (cognitive parameter)
- **c2 = 1.5**: Hệ số học tập xã hội (social parameter)
- **num_particles = 30**: Số lượng hạt trong bầy đàn
- **max_iterations = 100**: Số vòng lặp tối đa

## Cách chạy

### 1. Chạy thuật toán PSO:
```matlab
test_SIM_PSO
```

### 2. So sánh với thuật toán gốc:
```matlab
% Chạy thuật toán gradient descent gốc
test_SIM

% So sánh kết quả
load('Capacity_K_10.mat');      % Original
load('Capacity_K_10_PSO.mat');  % PSO

figure;
plot(1:10, Capacity_K_10, 'b-o', 'LineWidth', 2, 'DisplayName', 'Gradient Descent');
hold on;
plot(1:10, Capacity_K_10_PSO, 'r-s', 'LineWidth', 2, 'DisplayName', 'PSO');
xlabel('Number of TX-SIM Layers (L)');
ylabel('Sum Rate (bits/s/Hz)');
title('Comparison: Gradient Descent vs PSO');
legend('Location', 'best');
grid on;
```

## Đặc điểm của thuật toán PSO

### Ưu điểm:
1. **Global search**: PSO có khả năng tìm kiếm toàn cục tốt hơn gradient descent
2. **Không cần đạo hàm**: Không yêu cầu tính gradient của hàm mục tiêu
3. **Đơn giản**: Dễ hiểu và dễ implement
4. **Song song hóa**: Có thể đánh giá nhiều hạt song song

### Nhược điểm:
1. **Chi phí tính toán cao**: Cần đánh giá nhiều điểm trong không gian tìm kiếm
2. **Hội tụ chậm**: Có thể cần nhiều vòng lặp hơn gradient descent
3. **Điều chỉnh tham số**: Hiệu quả phụ thuộc vào việc chọn tham số w, c1, c2

## Kết quả mong đợi

File output:
- `NMSE_K_10_PSO.mat`: Normalized Mean Square Error theo số lớp L
- `Capacity_K_10_PSO.mat`: Sum Rate (Spectral Efficiency) theo số lớp L

Các biểu đồ:
1. NMSE vs Number of TX-SIM Layers
2. Spectral Efficiency vs Number of TX-SIM Layers

## Tham khảo

**Original Paper:**
J. An et al., "Stacked Intelligent Metasurfaces for Efficient Holographic MIMO Communications in 6G," 
IEEE Journal on Selected Areas in Communications, vol. 41, no. 8, pp. 2380-2396, Aug. 2023

**PSO Algorithm:**
Kennedy, J., & Eberhart, R. (1995). Particle swarm optimization. 
Proceedings of ICNN'95 - International Conference on Neural Networks.
