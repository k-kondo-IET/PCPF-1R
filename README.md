# PCPF family

Simulation models for predicting fate and transport of pesticide in paddy environment.

## Main functions

-   `PCPF-1R.R` : Plot scale simulation model consisting of paddy water and 1cm thick soil
-   `PCPF-BR.R` : Block scale simulation to calculate pesticide concentration in drainage canal

The codes include following sub-functions

-   `PPCPF.R` : Parameter estimation function using raw data
-   `PPONC.R` : Pre-processing observed concentration data based on FOCUS guidance
-   `WBcalc.R` : Function to calculate water balance from specified parameters
-   `FAO56ET.R` : Calculate daily evapotranspiration using FAO56 Penmann-Monteith method

## Requirements

-   [R](https://cran.r-project.org/)
-   [RStudio](https://www.rstudio.com/)

## Usage

1.  Select "Download ZIP" from "Code".
2.  Unzip the file to your local working directory.
3.  Open `PCPF-1R.Rproj`.
4.  Test `Hands-on.R` in "analysis" folder.

## Author

Kei Kondo, PhD

The Institute of Environmental Toxicology

Mail : kondoh[at]iet.or.jp

## References

1.  Kondo, K., Wakasone, Y., Iijima, K., Ohyama, K., 2019. Inverse analysis to estimate site-specific parameters of a mathematical model for simulating pesticide dissipations in paddy test systems. Pest. Manag. Sci., 75(6): 1594-1605. <https://doi.org/10.1002/ps.5276>
2.  Kondo, K., Wakasone, Y., Iijima, K., Ohyama, K., 2020. Inverse modeling of laboratory experiment to assess parameter transferability of pesticide environmental fate into outdoor experiments under paddy test systems. Pest Manag Sci, 76(8): 2768-2780. <https://doi.org/10.1002/ps.5824>
