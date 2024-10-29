# Geodynamic Activity Characteristics in the Sichuan–Yunnan Region Revealed by Big Data Analysis

## 0. Description
This document describes a new visualization methods for geoscience volume data, and utilizes the CVM of "SWChinaCVM-2.0" for the Sichuan–Yunnan region, combining scientific visualization methods such as Ray Casting volume rendering, Marching Cubes Iso-surface extraction, and hybrid rendering. We designed a two-dimensional transfer function related to velocity magnitude to quantitatively express high and low-velocity anomalies. By integrating geoscience big data such as high-resolution topography, three-dimensional active faults, earthquake catalogs, and the Moho discontinuity, we visually and quantitatively displayed the spatial distribution and thickness of velocity anomalies at different layers within the CVM. Additionally, we illustrated their relative relationships with faults, the Moho discontinuity, and historical large earthquake hypocenters using three-dimensional stereograms and animations with varying camera perspectives.

## 1. Background
The Sichuan–Yunnan region, located at the southeastern margin of the Tibetan Plateau,GPS observations indicate that present-day crustal movement in this region predominantly exhibits a clockwise rotational pattern.
<p align="center">
  <img src="pic/pic1.png" alt="Tectonic background map of the Sichuan–Yunnan region and its surroundings" /><br />
  Tectonic background map of the Sichuan–Yunnan region and its surroundings
</p>

It's characterized by complex tectonic deformation, deep active faults, and frequent earthquakes.
<p align="center">
  <img src="pic/pic2.png" alt="2D/3D scatter plot of earthquakes, 3D faults and elevation" /><br />
  2D/3D scatter plot of earthquakes, 3D faults and elevation
</p>
So, the establishment of community velocity model(CVM) is as important as their visualization and post-processing analysis in this region.

## 2. Data
#### (1) [SWChinaCVM-2.0 Model Data](http://cses.ac.cn/sjcp/ggmx/2022/589.shtml);
#### (2) [High-Resolution Topography Data (NOAA National Centers for Environmental Information, 2022)](https://www.ncei.noaa.gov/products/etopo-global-relief-model);
#### (3) [Three-Dimensional Fault Data](http://cses.ac.cn/sjcp/ggmx/2024/609.shtml);
#### (4) [Crustal Depth Data](https://www.sciencedirect.com/science/article/abs/pii/S0040195113006847?via%3Dihub);
#### (5) [GPS Velocity Field](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019JB018774);
#### (6) Earthquake Catalog Data.

## 3. Method
### 3.1 Transfer Functions of Velocity
In geosciences, the velocity anomalies of interest typically refer to regions that are significantly higher (high-velocity anomalies) or lower (low-velocity anomalies) than the background velocity values within the same depth range. Since velocity generally increases with depth, the mapping method of scalar velocity values to colors varies for each layer. Therefore, it is necessary to unify the mapping range of the entire CVM volume data first:
<p align="center">
  <img src="pic/tf_func.png" alt="Transfer Function of Velocity" /><br />
</p>

### 3.2 Comparison
<p align="center">
  <img src="pic/pic3.png" alt="" /><br />
  Comparison of the separation effect between the CVM data and high- and low-velocity anomalies mapped by the Transfer Functions.
</p>

<p align="center">
  <img src="pic/pic4.png" alt="" /><br />
  Comparison between the 3D reconstruction results of velocity anomalies and the original CVM data. 
</p>

## 4. Results
### 4.1 3D spatial distribution of the Vp anomalies in the Sichuan–Yunnan region
<p align="center">
  <img src="pic/vp.png" alt="" /><br />
  3D spatial distribution of P wave velocity anomalies in the Sichuan–Yunnan region (vertical exaggeration of ten times). (a) High-velocity anomalies; (b) Low-velocity anomalies. Solid lines represent active faults. The light gray surface denotes the Moho surface.
</p>

### 4.2 3D spatial distribution of the Vs anomalies in the Sichuan–Yunnan region
<p align="center">
  <img src="pic/vs.png" alt="" /><br />
  3D spatial distribution of S wave velocity anomalies in the Sichuan–Yunnan region (vertical exaggeration of ten times). (a) High-velocity anomalies; (b) Low-velocity anomalies. Solid lines represent active faults. The light gray surface denotes the Moho surface.
</p>

### 4.3 3D spatial distribution of the Vp/Vs ratio anomalies in the Sichuan–Yunnan region
<p align="center">
  <img src="pic/vpDvs.png" alt="" /><br />
  Anomalous distribution of the Vp/Vs ratio in the Sichuan–Yunnan region (vertical exaggeration of ten times). (a) High V_p/V_s anomalies; (b) Low V_p/V_s anomalies. Solid lines represent active faults. The light gray surface denotes the Moho surface.
</p>

### 4.4  Volume Rendering of the Vp and Vs Anomalies
<p align="center">
  <img src="pic/pic4.png" alt="" /><br />
  Volume rendering results of the Vp(left) and Vs(right) combining high- and low-velocity anomalies (vertical exaggeration of ten times).Warm red indicates low-velocity anomalies; cool green indicates high-velocity anomalies. Solid lines represent active faults.
</p>

### 4.5  Volume Rendering of the Vp/Vs Ratio Anomalies
<p align="center">
  <img src="pic/pic4.png" alt="" /><br />
  Volume rendering results of the Vp/Vs ratio combining high- and low-velocity anomalies (vertical exaggeration of ten times).Warm red indicates low-velocity anomalies; cool green indicates high-velocity anomalies. Solid lines represent active faults. 
</p>

### 4.6 Comprehensive Spatial Distribution Relationship Diagram
<p align="center">
  <img src="pic/vp_vs_10km.png" alt="" /><br />
  Spatial relationships between high- and low-velocity anomalies and earthquake hypocenter locations at different depth ranges (vertical exaggeration of ten times). In each panel, the left side shows P wave anomalies, and the right side shows S wave anomalies. Warm red indicates low-velocity anomalies; cool green indicates high-velocity anomalies. Colored vertical planes represent 3D active faults. Small particle spheres represent earthquakes within the respective depth ranges. The light gray surface represents the Moho discontinuity.
</p>

## 5. Comprehensive 3D Animation
This video integrates over 420,000 earthquake entries of magnitudes 1.0 and above from the Sichuan-Yunnan region, covering the period from 1970 to 2020. It includes data on more than 160 three-dimensional active faults, high-precision topographic information, community velocity models, and the depth of the Moho surface, which comprehensives analysis reveals the geodynamic activity characteristics of the Sichuan–Yunnan region.

[![Watch the video](https://img.youtube.com/vi/x29ss0pleRU/maxresdefault.jpg)](https://youtu.be/x29ss0pleRU)
*(Complete 3D Visualization Video: [Visit on YouTube](https://youtu.be/x29ss0pleRU))*



## How to Use
### 1. Installation

### 1.1 Clone the Repository
Clone the repository to your local machine:
```bash
git clone https://github.com/peixyhe/seismic-catalog_vis_demo.git
cd seismic-catalog_vis_demo-master
```

### 1.2 Create a Conda Environment
This project relies on several libraries, which can be installed via the `environment.yml` file. Key dependencies include:
- **VTK**: For 3D visualization and rendering.
- **Pandas**: For data manipulation and analysis.
- **NumPy**: For numerical operations.
- **SciPy**: For spatial data processing with `cKDTree`.
- **tqdm**: For progress bars.
- **joblib**: For parallel computation.

Install all required dependencies:
```bash
conda env create -f environment.yml
conda activate seismic_vis_env
```

### 1.3 Usage
Run the desired script. If there are any unclear areas, please refer to the comments regarding the main function in each script. For example:

#### Example 1: 3D Earthquake Catalog Scatter Plot
```bash
python ./src/catalog_scatter_vis.py 3D ./rawData/CENC_catalog_1970-2020.csv
```
This command inputs an earthquake catalog CSV file and outputs a PNG result and a VTK file for further visualization.

#### Example 2: M-T Heat Map
```bash
python ./src/MT_hotmap_vis.py ./rawData/CENC_catalog_1970-2020.csv
```
This command also inputs an earthquake catalog CSV file, generating a PNG result and a VTK file.

## Citation
If this work is useful to you, please cite the following source: He Pei et al., "Geodynamic Activity Characteristics in the Sichuan–Yunnan Region Revealed by Big Data Analysis" This work is currently under review, and specific publication details will be provided upon acceptance.

## License
This project is licensed under the MIT License. See the [MIT License](LICENSE) file for details.

---
