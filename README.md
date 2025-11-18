### Project Description

GPS Satellites Constellation is a visualization and simulation project that shows how GPS satellites move around Earth and how they help determine your location.

The project displays:
- The 3D orbits of GPS satellites in space  
- Their real-time positions relative to Earth  
- Which satellites are visible from any chosen location on the ground (using a 7° elevation mask)  
- How the number of visible satellites changes over time  
- Simulated signal distances (pseudoranges), including:
  - Real geometric distance  
  - A fixed clock bias (+1.5 m)  
  - Random measurement noise (±3 m)  
- The basic idea of trilateration - how GPS finds your position using signals from multiple satellites
- A positioning algorithm was implemented to calculate the receiver’s coordinates from the simulated pseudoranges. The resulting position error over time is shown in the final graph.

This tool is designed for students, hobbyists, and anyone who wants to understand how GPS actually works — without complex math or real hardware.

---

### Result Visualization

#### 1. 3D Orbital Constellation  
<img width="600" alt="3D orbits" src="https://github.com/user-attachments/assets/59a71e22-4f5a-4de9-9c32-db1476f87983" />

#### 2. Number of Visible Satellites Over Time  
<img width="600" alt="Visible satellites" src="https://github.com/user-attachments/assets/0ce3c532-08e3-4972-9125-e9494a60f318" />

#### 3. Simulated Pseudoranges for Five Satellites  
<img width="600" alt="Pseudoranges" src="https://github.com/user-attachments/assets/d068a78e-59d0-4e99-b98b-2b23ed9e55e9" />

#### 4. Positioning Error Over Time  
<img width="600" alt="изображение" src="https://github.com/user-attachments/assets/e00af9e2-de6a-42df-953f-3e75e391a616" />


