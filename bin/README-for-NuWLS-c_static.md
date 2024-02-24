## NuWLS-c_static README

### Introduction
NuWLS-c_static is a solver for (W)PMS. This README provides instructions on using NuWLS-c_static, including optional parameters and examples.

### Usage
NuWLS-c_static supports an optional parameter `-init_soln_file filename`. This parameter specifies a file (`filename`) containing the initial assignment for the SLS solver NuWLS. If this parameter is provided, the initial assignment from the file will be used.

### File Format
The `filename` should follow the format of the file `init_assignment_ram_k3_n20.ra1.txt`.

### Example Run
```bash
./NuWLS-c_static ram_k3_n20.ra1.wcnf -init_soln_file init_assignment_ram_k3_n20.ra1.txt
```
or
```bash
./NuWLS-c_static ram_k3_n20.ra1.wcnf
```

### Running with Scripts
Alternatively, use scripts for execution. For example:
```bash
./starexec_run_default-runsolver-NuWLS2023 ram_k3_n20.ra1.wcnf -init_soln_file init_assignment_ram_k3_n20.ra1.txt
```
or
```bash
./starexec_run_default-runsolver-NuWLS2023 ram_k3_n20.ra1.wcnf
```
or
```bash
./starexec_run_short-runsolver-NuWLS2023 ram_k3_n20.ra1.wcnf -init_soln_file init_assignment_ram_k3_n20.ra1.txt
```
or
```bash
./starexec_run_short-runsolver-NuWLS2023 ram_k3_n20.ra1.wcnf
```

This README provides basic usage instructions for NuWLS-c_static. For more detailed information, please refer to the documentation or contact the maintainers.
