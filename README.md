# Replication Code for CGLS (2022)

This is the replication code for "Do Low Interest Rates Harm Innovation, Competititon, and Productivity Growth," by Craig A. Chikis, Jonathan Goldberg, and David López-Salido (2022). It contains codes to replicate our numerical results. 

*Abstract*: Real interest rates, productivity growth, and business dynamism have fallen in recent decades.  To assess whether these trends are driven by low financial discounts, we develop and estimate a Schumpeterian model that nests alternative assumptions about creative destruction and includes strategic behavior among innovating firms.  Our framework provides a good fit to many cross-sectional patterns related to firm growth, innovation, and markups.  We find that low discount rates increase productivity growth and, to some extent, market competition.  Low discount rates can entrench low-R&D market leaders if creative destruction is severely restricted by assumption.  However, such severe restrictions imply extremely and counterfactually weak business dynamism.  Many additional exercises, including restricting laggards' innovativeness in uncompetitive industries, illustrate the robustness of our results.  Our findings suggest that recent productivity and dynamism trends are not explained by low financial discounts.  

 
## Getting Started

The first thing to do is to clone this repository. Remember where you save it, because you will need to tell the main script where to find it in your local file system. 

```	
git clone https://github.com/cchikis/CGLS_RepCode.git
```

If you'd rather not work via the command line, you can simply download a `.zip` file with the codes. 

### External Dependencies

We use NLopt to generate Figure IA.16.  If NLopt is not installed, the code runs successfully but does not generate this figure.  Instructions on how to install NLopt can be found [here](https://nlopt.readthedocs.io/en/latest/).

You will also need the MATLAB Parallel Computing and Optimization Toolboxes. Furthermore, it will be helpful to have a machine with access to multiple CPUs.

This replication package contains codes that replicate our numerical results. The code that generates our empirical results is also included in `Code/Data`, but the data itself is not, because some of the data is subject to data use agreements that prohibit redistribution.  All of the data used is either publicly available (see Internet Appendix for details and links) or available in commonly accessible databases (e.g., WRDS/Compustat).

## Generating All Figures and Tables

Go to `Code/Main/main_cgls.m`. As a string or character-type, enter the directory where you saved the repository for `main_dir`, and the directory where you saved the MATLAB API for NLopt for `nlopt_dir`. For me, operating on a Linux system, it looks like this.

``` MATLAB
% Fill this in per your local file setup
main_dir = "/home/cachikis/CGLS_RepCode/"; 

% You need NLopt for some of this code; enter the path where the MATLAB API is saved
% https://nlopt.readthedocs.io/en/latest/ Has instructions how to install (it is open-source)
nlopt_dir = "/home/cachikis/NLopt/mex/";
% If you don't have it just uncomment this line
% nlopt_dir = ""; 
```

Once you've done that, you can simply run `Code/Main/main_cgls.m`. The function at the bottom, `run_paper_new();` will execute the necessary files. Figures with appear in `Output/Figures_Paper/` and tables will appear in `Output/LaTeX_Output/`. 

