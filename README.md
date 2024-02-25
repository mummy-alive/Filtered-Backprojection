# Filtered-Backprojection
2023-11-20 Medisys lab intern study

### Description
This is a Matlab code of Filtered Backprojection of Circle $$(x-a)^2 + (y-b)^2 = r^2$$ on coordinate plane.\

Input
- Coordinate value of center of the circle x, y
- Radius
- Number of detectors & Sources (This shall an odd number)
- Distance between each detectors
- SDD, SCD
- Total views of CT
- Attenuation coefficient of circle

### Step 1: Forward Projection
testcase 1:
```
5
3
6
101
0.5
20
9
30
1
```
testcase 2:
```
5
3
12
101
0.5
30
14
72
1
```
output:![Untitled](https://github.com/mummy-alive/Filtered-Backprojection/assets/113423544/3c206e72-870b-4b17-852c-6b19a51c2567)

### Step 2. Fourier Transform & Filter
output: ![Filtered sinogram](https://github.com/mummy-alive/Filtered-Backprojection/assets/113423544/30a5ffc8-631c-41e5-a79b-9e13a84926cf)
Note that the Filtered sinogram graph requires resizing.

### Step 3. Interpolation & Inverse Fourier Transform & FBP

Result: ![FBP result](https://github.com/mummy-alive/Filtered-Backprojection/assets/113423544/4c832f2f-4b2e-48ac-9f32-c30b8053e051)
