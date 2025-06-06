# LidDrivenCavity
This is an implementation of CFD codes for the Lid Driven Cavity using OpenMP, C++ and Vorticity Streamfunction formulation.

* [Paper](https://www.iist.ac.in/sites/default/files/people/psi-omega.pdf)
* [Video](https://www.youtube.com/watch?v=3JnBTB0Edak)
* [Report](https://github.com/user-attachments/files/20624129/Tau_Class___Lab_Report_Template.1.pdf)



## To run
- Step 0;
  ```bash
  git clone https://github.com/varun-r-mallya/LidDrivenCavity.git;
  cd LidDrivenCavity; git submodule update --init --recursive`
  ```
- Step 1: `meson setup buildDir`
- Step 2: `cd buildDir`
- Step 3: `ninja`
- Step 4: `./LidDrivenCavity /path/to/your/config/file`

### Configuration
```toml
[LidDrivenCavity]
Re = 500.0              # Reynolds number
length = 2.0            # Length of Grid
grid_points = 100       # Number Of GridPoints
time_step = 0.01        # Time step
max_iterations = 10     # Maximum Number of iterations
tolerance = 1e-6        # tolerance adjust
python_plot = true
```


