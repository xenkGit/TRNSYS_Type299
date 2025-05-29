# TRNSYS_Type299
TRNSYS model for irregular borehole heat exchanger fields.

## Description
This model represents a double U-tube borehole heat exchanger field with irregular geometry, coupling a short-term thermal resistance and capacitance model and a combination of g-functions. The coordinates of the boreholes should be specified in a file called "coordinates.txt" and placed in the same folder as the dck-file.

By default, the g-function is calculated internally. Alternatively, a g-function can be provided in a file called "pygfunction_gFunc.txt" (e.g., generated with pygfunction [1]) and placed in the same folder as the dck-file (in this case make sure to use the same time settings in pygfunction and cells_per_level=5).

The vertical and horizontal discretisation can be set using "Number of segments" and "Number of sheath grout nodes". Please refer to the documentation for the model limitations.

## Dependencies
The DLLs were compiled for TRNSYS 18.06 (64 bit) using the v142 toolset and C++20 standard.

To compile the C++ source code, you have to include and link the GNU scientific library (GSL) v2.3.0.2779 [2] and include the Eigen library v3.4.0 [3].

## References
[1] https://pygfunction.readthedocs.io/en/stable/# <br>
[2] https://www.gnu.org/software/gsl/ <br>
[3] https://eigen.tuxfamily.org/index.php?title=Main_Page
