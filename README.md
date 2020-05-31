# KWAVE

This model is designed to represent infiltration (Green-Ampt), interception (Rutter), and runoff (kinematic wave).
All the model code is contained in "kinematic.c"
Matlab scripts are used to write model input files and read model output files

The easiest way to run the example code, which would require Matlab, is to:
(1) Execute "WriteInputFiles.m" to write the input files that the model will be looking for
(2) Compile and Execute "kinematic.c" to run the model
(3) When the C code had finished running, execute "ReadOutputFiles.m" to plot results 

If you follow these steps to run the example as it is set up, the C code should finish 
running in several minutes. The example code is set up to simulate runoff in a small
watershed in response to a short (15-minute duration) rainstorm. 

Alternatively, you can construct the input files in whatever way is convenient for you.
You will need to make sure that they have the same format as those that are produced by
the Matlab code "WriteInputFiles.m".


For details on the model and numerical solution, see the following references:

Rengers, F.K., McGuire, L.A., Kean, J.W., Staley, D.M. and Hobley, D.E.J., 2016. 
Model simulations of flood and debris flow timing in steep catchments after wildfire. 
Water Resources Research, 52(8), pp.6041-6061.

McGuire, L.A. and Youberg, A.M., 2019. Impacts of successive wildfire on soil hydraulic properties: 
Implications for debris flow hazards and system resilience. Earth Surface Processes and Landforms, 
44(11), pp.2236-2250.

Input Files:
1. input -- must 1 row with tab seperated values for the following parameters arranged in order: 
[nx ny dx tend epsilon h0 rnum1 rint d84 pi Si Ki gi]
See the matlab script for parameter explanations
2. topoin -- a grid of elevation data, reshaped into a single row with tab seperated values 
3. depthin -- a grid of initial flow depth, reshaped into a single row with tab seperated values 
4. solidin -- a grid containing 0s in all locations within the coputational domain and 1 in areas outside of the computational domain
5. rain1 -- times series of rainfall data
5. rain2 -- times series of rainfall data
7. ksin -- a grid of saturated hydraulic conductivity, reshaped into a single row with tab seperated values 
8. vinfin -- a grid with the depth of water infiltrated prior to start of simulation, reshaped into a single row with tab seperated values 
9. vegcoverin -- a grid of vegetation cover fraction, reshaped into a single row with tab seperated values 
10. channelin -- a grid with a 1 in all locations identified as a channel and 0 otherwise, reshaped into a single row with tab seperated values 

Output Files:
1. stage -- array containing information on flow at the edges of the model domain
2. depth -- flow depth at each grid cell at the end of the simulation
3. vel -- flow velocity at each grid cell at the end of the simulation
4. maxdepth -- maximum flow depth at each grid cell
4. maxvel -- maximum flow velocity at each grid cell

