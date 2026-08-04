---
name: geoskill-precipitation-nowcasting
description: '基于光流法（交叉相关位移估计）的拉格朗日持久性降水临近预报，外推未来 0-6 小时降水场，输出预报序列 GeoTIFF 与位移场 JSON。Optical-flow (cross-correlation) Lagrangian persistence nowcasting that extrapolates precipitation fields 0-6 hours ahead, outputting a forecast GeoTIFF stack and a displacement-field JSON.'
---

# 降水临近预报 | Precipitation Nowcasting

An optical-flow based Lagrangian persistence precipitation nowcast that estimates the motion of the precipitation field from the most recent radar / satellite frames and extrapolates it 0–6 hours ahead. Suitable for very-short-range precipitation forecasting, heavy-rain warnings, and pre-assessment of urban flooding and flash floods.

Core algorithm:

- **Displacement estimation (optical flow)**: matches pairs of adjacent precipitation fields within a search window using normalized cross-correlation to estimate the translation vector (vy, vx); results are averaged over multiple adjacent frame pairs for robustness (simplified cross-correlation optical flow, equivalent to Lucas-Kanade under the assumption of uniform translation).
- **Lagrangian extrapolation**: assuming the displacement field is stationary over short lead times, the latest frame is shifted along the estimated velocity (bilinear resampling via scipy.ndimage.shift, zero-padding at the borders) to produce the forecast fields for each future lead time.

A built-in `--synthetic` mode generates a sequence of Gaussian rain cells translated at known velocities, for offline validation of the displacement estimation and extrapolation accuracy.

## Dependencies / 依赖

```bash
pip install numpy rasterio scipy
```

## Usage / 使用方法

### Basic usage (synthetic data, offline)

```bash
python geoskill-precipitation-nowcasting.py --bbox 116.0 39.0 117.0 40.0 --output-dir ./output
```

### Example 1: forecast 60 minutes ahead (4-frame input)

```bash
python geoskill-precipitation-nowcasting.py --bbox 116 39 117 40 --n-frames 4 --lead-time 60 --output-dir ./nowcast_60
```

### Example 2: forecast 90 minutes ahead, 10-minute steps

```bash
python geoskill-precipitation-nowcasting.py --bbox 121 31 122 32 --lead-time 90 --dt-minutes 10 --output-dir ./nowcast_90
```

### Example 3: real radar sequence (multi-band, one band per frame)

```bash
python geoskill-precipitation-nowcasting.py --input radar_stack.tif --lead-time 60 --search 16 --output-dir ./radar_nowcast
```

### Example 4: larger search radius (fast-moving systems)

```bash
python geoskill-precipitation-nowcasting.py --bbox 116 39 117 40 --n-frames 5 --search 20 --lead-time 120 --output-dir ./fast_system
```

### Example 5: bbox-only auto synthesis + quiet mode

```bash
python geoskill-precipitation-nowcasting.py --bbox 110 30 111 31 --lead-time 60 --output-dir ./auto --quiet
```

## Output / 输出

| File | Format | Description |
|---|---|---|
| `forecast.tif` | GeoTIFF (float32, N band) | Extrapolated precipitation field per lead time; band count = number of forecast steps, EPSG:4326 |
| `displacement.json` | JSON | Mean displacement vector, per-frame-pair displacements, peak correlation, forecast lead times |
| `output-manifest.json` | JSON | Run manifest (inputs/outputs/QA/software versions) |

## Data Source / 数据源 / Source

- **Input mode**: local multi-band GeoTIFF (one band per time step of the precipitation field).
- **Synthetic mode**: locally generated translated Gaussian rain cells; no external data sources.

## Privacy / 隐私声明 / Privacy

- Fully offline by default; no network requests are made.
- `--synthetic` mode reads no external data.
- All computation is performed locally; user data is never uploaded.

## License / License

MIT

---

<!-- ===== 中文原文 (Chinese Original) ===== -->

---
name: geoskill-precipitation-nowcasting
description: '基于光流法（交叉相关位移估计）的拉格朗日持久性降水临近预报，外推未来 0-6 小时降水场，输出预报序列 GeoTIFF 与位移场 JSON。Optical-flow (cross-correlation) Lagrangian persistence nowcasting that extrapolates precipitation fields 0-6 hours ahead, outputting a forecast GeoTIFF stack and a displacement-field JSON.'
---

# 降水临近预报 | Precipitation Nowcasting

基于**光流法**的拉格朗日持久性（Lagrangian persistence）降水临近预报，
利用最近几帧雷达 / 卫星降水场估计场体移动速度，并向未来 0–6 小时外推。
适用于短临降水预报、暴雨预警、城市内涝与山洪的前置研判。

核心算法：

- **位移估计（光流）**：在搜索窗内用归一化互相关匹配相邻两帧降水场，
  估计平移矢量 (vy, vx)；对多个相邻帧对取平均以提高稳健性（简化交叉相关
  光流，等价于全场平移假设的 Lucas-Kanade）。
- **拉格朗日外推**：假设位移场短时不变，将最新一帧沿估计速度平移
  （scipy.ndimage.shift 双线性重采样，边界外补零），得到未来各时次预报场。

内置 `--synthetic` 模式生成以已知速度平移的高斯雨团序列，用于离线验证
位移估计与外推精度。

## 依赖

```bash
pip install numpy rasterio scipy
```

## 使用方法

### 基本用法（合成数据，离线）

```bash
python geoskill-precipitation-nowcasting.py --bbox 116.0 39.0 117.0 40.0 --output-dir ./output
```

### 示例 1：预报未来 60 分钟（4 帧输入）

```bash
python geoskill-precipitation-nowcasting.py --bbox 116 39 117 40 --n-frames 4 --lead-time 60 --output-dir ./nowcast_60
```

### 示例 2：预报未来 90 分钟，10 分钟步长

```bash
python geoskill-precipitation-nowcasting.py --bbox 121 31 122 32 --lead-time 90 --dt-minutes 10 --output-dir ./nowcast_90
```

### 示例 3：真实雷达序列（多波段，每波段 = 一帧）

```bash
python geoskill-precipitation-nowcasting.py --input radar_stack.tif --lead-time 60 --search 16 --output-dir ./radar_nowcast
```

### 示例 4：更大搜索半径（快速移动系统）

```bash
python geoskill-precipitation-nowcasting.py --bbox 116 39 117 40 --n-frames 5 --search 20 --lead-time 120 --output-dir ./fast_system
```

### 示例 5：仅 bbox 自动合成 + 静默

```bash
python geoskill-precipitation-nowcasting.py --bbox 110 30 111 31 --lead-time 60 --output-dir ./auto --quiet
```

## 输出

| 文件 | 格式 | 说明 |
|---|---|---|
| `forecast.tif` | GeoTIFF (float32, N band) | 逐时次外推降水场，band 数 = 预报步数，EPSG:4326 |
| `displacement.json` | JSON | 平均位移矢量、逐帧对位移、峰值相关、预报时次 |
| `output-manifest.json` | JSON | 运行清单（输入/输出/QA/软件版本） |

## 数据源 / Source

- **输入模式**：本地多波段 GeoTIFF（每波段 = 一个时间步的降水场）。
- **合成模式**：本地生成平移高斯雨团，无外部数据源。

## 隐私声明 / Privacy

- 默认完全离线运行，不发起任何网络请求。
- `--synthetic` 模式不读取任何外部数据。
- 所有计算在本地完成，不上传用户数据。

## License

MIT
