Вот обновлённое описание с **единым, средним размером изображений** (ширина `600px` — достаточно для читаемости, но не занимает всю ширину экрана):

---

### Project Description

**GPS Satellites Constellation** is a visualization and simulation project that demonstrates the orbital mechanics and positioning of satellites in the Global Positioning System (GPS) constellation. The project aims to provide an intuitive understanding of how GPS satellites orbit the Earth, their distribution in space, and how they enable accurate positioning and navigation on the planet's surface.

Using simplified models and visualizations, this project shows:
- The orbital paths of GPS satellites in 3D inertial space  
- Satellite positions at any given time, transformed into Earth-fixed coordinates (DGSK)  
- Visibility analysis from a user-defined ground location (latitude, longitude, altitude)  
- **Time-varying number of visible satellites**, accounting for elevation mask angle (7°)  
- **Simulated pseudoranges** for visible satellites, incorporating:
  - Geometric distance from receiver to satellite  
  - Constant receiver clock bias (**+1.5 m**, as per lab guidelines)  
  - Uniform measurement noise (**±3 m**)  
- Basic principles of trilateration used in GPS calculations  

This project is ideal for educational purposes, space enthusiasts, and developers interested in satellite technology and orbital mechanics.

---

### Result Visualization

#### 1. 3D Orbital Constellation  
<img width="600" alt="3D orbits" src="https://github.com/user-attachments/assets/59a71e22-4f5a-4de9-9c32-db1476f87983" />

#### 2. Number of Visible Satellites Over Time  
<img width="600" alt="изображение" src="https://github.com/user-attachments/assets/0ce3c532-08e3-4972-9125-e9494a60f318" />


#### 3. Simulated Pseudoranges for Five Satellites  
<img width="600" alt="изображение" src="https://github.com/user-attachments/assets/d068a78e-59d0-4e99-b98b-2b23ed9e55e9" />
