# LPT Calibration System — README
## 4-Camera Refractive Lagrangian Particle Tracking Calibration (MATLAB)

---

## Overview

This package implements a complete calibration pipeline for a 4-camera
3-D Lagrangian Particle Tracking (LPT) system in a water flume, correctly
accounting for refraction at the **air → glass → water** interfaces.

---

## File Structure

```
LPT_Calibration/
│
├── LPT_Config.m                ← START HERE — set all system parameters
│
├── Phase1_AirCalibration.m     ← Intrinsics (focal length, distortion) via checkerboard in air
├── Phase2_WandDetection.m      ← LED wand detection + interactive verification GUI
├── Phase3_RefractionModel.m    ← Core refractive ray model (utility, not run directly)
├── Phase4_BundleAdjustment.m   ← Nonlinear bundle adjustment for extrinsic poses
├── Phase5_Reconstruction.m     ← Refraction-corrected 3-D particle triangulation
│
├── Visualize_System.m          ← 3-D visualisation of camera network + ray paths
├── EvaluateCalibration.m       ← Quantitative accuracy assessment
├── ExportCalibration.m         ← Export to JSON / CSV / OpenPTV format
├── SyntheticDataGenerator.m    ← Generate synthetic data to validate pipeline
│
└── data/
    ├── air_calib/
    │   ├── cam1/   ← checkerboard images for camera 1 (in air)
    │   ├── cam2/
    │   ├── cam3/
    │   └── cam4/
    └── wand_images/
        ├── cam1/   ← synchronised wand frames for camera 1 (in water)
        ├── cam2/
        ├── cam3/
        └── cam4/
```

---

## Recommended Workflow

### Step 0 — Validate with synthetic data (do this first)
```matlab
SyntheticDataGenerator    % generate synthetic wand observations
Phase4_BundleAdjustment   % calibrate using synthetic data
EvaluateCalibration       % compare to ground truth
```
If the synthetic test fails, there is a bug in your config or code before
you waste time collecting real data.

### Step 1 — Configure
Edit `LPT_Config.m`:
- Set `imageSize` for your Edgertronic model
- Set `walls(c).normal`, `.point`, `.thickness` for each camera's glass wall
- Set `wandLength` (measure carefully with calipers!)
- Set `airCalibDir` and `wandImageDir` to your data folders

### Step 2 — Air calibration (intrinsics)
```matlab
Phase1_AirCalibration
```
Collect 40–80 checkerboard images per camera **in air** (flume empty or cameras
pointed at a dry surface). Cover all regions of the image at multiple distances
and tilt angles. Target mean reprojection error < 0.3 px.

### Step 3 — Wand detection
```matlab
Phase2_WandDetection
```
After filling the flume and setting up illumination, wave the LED wand
through the full measurement volume. Aim for ≥500 frames where all 4
cameras see the wand. Run this script, then use the GUI to reject bad
detections (arrow keys to navigate, 'T' to toggle frame validity).

### Step 4 — Bundle adjustment
```matlab
Phase4_BundleAdjustment
```
This is the core calibration step. Runtime depends on number of frames
(~5–20 minutes for 500 frames on a modern CPU). Target final RMS
reprojection error < 0.5 px.

### Step 5 — Evaluate & visualise
```matlab
EvaluateCalibration
Visualize_System
```

### Step 6 — Export
```matlab
ExportCalibration
```

### Step 7 — Reconstruct particles
```matlab
Phase5_Reconstruction
```
Provide your particle observations in `results/particle_observations.mat`
(format described inside the script).

---

## Physical Model

Each camera ray passes through **two planar refractive interfaces**:

```
Camera (air)
    │
    │  d_air
    ▼
╔═════════════╗  ← outer glass face
║    GLASS    ║  ← d_glass (Snell: n_air→n_glass)
╚═════════════╝  ← inner glass face
    │
    │  d_water (Snell: n_glass→n_water)
    ▼
  Water (measurement volume)
    ·  ← particle
```

**Snell's law (vector form):**

```
d_out = (n1/n2)*d_in + [(n1/n2)*cos(θ_i) - cos(θ_t)] * n̂
```

where `cos(θ_t) = sqrt(1 - (n1/n2)² * sin²(θ_i))`.

The wall geometry is defined per-camera by:
- `normal` — outward unit normal (pointing into air)
- `point`  — any point on the **inner** (water-side) face
- `thickness` — glass wall thickness

Cameras at **arbitrary angles** to the wall are fully supported. The
only geometric assumption is that the glass wall is **planar** (flat glass).

---

## LED Wand Design

### Recommended specification

| Parameter | Value |
|---|---|
| LED type | High-brightness IR (850 nm) or white |
| LED current | 50–100 mA (use current-limiting resistors) |
| Separation | 100 mm (measure with calipers, ±0.1 mm) |
| Wand material | Carbon fibre or aluminium tube (rigid) |
| LED mounting | Recessed cups to reduce glare spread |
| Encapsulation | Waterproof epoxy or heat-shrink tubing |

### Why IR LEDs?
- Invisible to the human eye → won't affect flow if particles are fluorescent
- Edgertronic cameras have good sensitivity in near-IR
- Narrow spectrum → easy to filter out background illumination with a bandpass filter

### Wiring
```
Battery (3.7V LiPo)  →  100Ω resistor  →  LED1  →  LED2  →  GND
                                          (series)
```
Mount the battery and switch in a waterproof housing at the handle end.

### Calibration tips
- Wave slowly (< 0.5 m/s) to minimise motion blur at typical frame rates
- Cover all corners and depths of the measurement volume
- Include wand orientations in all directions (don't just wave horizontally)
- Aim for ≥500 valid synchronised frames across all 4 cameras

---

## MATLAB Requirements

- MATLAB R2020b or later
- Computer Vision Toolbox (for `detectCheckerboardPoints`, `estimateCameraParameters`,
  `estimateEssentialMatrix`, `cameraIntrinsics`)
- Optimization Toolbox (for `lsqnonlin`)
- Image Processing Toolbox (for `bwconncomp`, `regionprops`, `imfilter`)

---

## Refractive Index Reference

| Medium | n (550 nm, 20°C) |
|---|---|
| Air | 1.0003 |
| Fresh water | 1.333 |
| Salt water (35 ppt) | 1.339 |
| Borosilicate glass | 1.474 |
| Soda-lime glass | 1.520 |
| Acrylic (PMMA) | 1.491 |
| Fused silica | 1.458 |

For salt water flumes, update `cfg.n_water` accordingly.
Temperature corrections: dn/dT ≈ −1.0×10⁻⁴ K⁻¹ for water.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| BA doesn't converge | Poor initialisation | Check that ≥3 cameras share frames; try more wand data |
| Reprojection error > 1 px | Wrong wall geometry | Double-check `wall.normal` and `wall.point` in Config |
| Scale is wrong | Wrong `wandLength` | Re-measure wand; check units (metres) |
| Cameras not found in images | Poor lighting | Use brighter LEDs; increase `intensityPct` threshold |
| TIR warnings in refraction | Ray at steep angle | Camera is looking at >~50° from wall normal — this is OK, rays at shallower angles still work |
| Phase 1 error > 0.5 px | Too few/poor images | Collect more checkerboard images at varied angles |

---

## Citation / Acknowledgement

If you use this code in a publication, please cite it as:

```
LPT Calibration System for Multi-Camera Water Flume Experiments
with Refractive Ray Model. MATLAB implementation. (2024)
```

---

*Contact: edit this line with your lab/project contact*
