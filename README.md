# GeoCalVal model
Calibration and Validation of a **linear model** created by Suhyb Salama, Department of Water Resources, ITC Faculty, University of Twente, The Netherlands
Based on *Salama, M. S., van der Velde, R.; Van der Woerd, H.J.; et al. (2012): Technical Notes: Calibration and validation of geophysical observation models. Biogeosciences 9,6: 2195-2201, DOI: 10.5194/bg-9-2195-2012*

# Way of operation:
## Input
Supply an input a file with a pair of X and Y;
The first column is the independent variables X; e.g. concentration/soil moisture. The second column is the dependent variables Y; e.g. absorption/ backscattering;
Example of a linear system is Y=slope * X+ intercept;
The idea is to find "slope" and "intercept" using the calibration set and estimate their accuracies using the validation set.
## Output:
Probability distribution functions of module coefficients (resulting from the calibration) 
Probability distribution functions of the error metrics, resulting from the validation 
# Note
It is straightforward to update for nonlinear model 
