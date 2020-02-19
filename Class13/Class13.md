Untitled
================

``` r
library(bio3d)
hsg_pdb<-get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

``` r
Hiv_prot<- read.pdb("1hsg.pdb")

Hiv_prot
```

    ## 
    ##  Call:  read.pdb(file = "1hsg.pdb")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

\#\#Q1: What is the name of the two non protein resid values in this
structure? What does resid correspond to and how would you get a listing
of all reside values in this structure? A: HOH and MK1 and residue is
the amino acids and you can get these values by listing or creating a
vector of the atoms-Resid list

\#\#Make ligand only ad prot only PDB Files

``` r
prot <- atom.select(Hiv_prot, "protein", value=TRUE)
lig <- atom.select(Hiv_prot, "ligand", value=TRUE)
write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, file="1hsg_ligand.pdb")
```

## Loading in and inspecting Docking Results

``` r
library(bio3d)
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```

``` r
res <- read.pdb("all.pdbqt", multi=TRUE)
ori <- read.pdb("ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1]  0.649  4.206 11.110 10.529  4.840 10.932 10.993  3.655 10.996 11.222
    ## [11] 10.567 10.372 11.019 11.338  8.390  9.063  8.254  8.978

``` r
rmsd(atom.select(ori,"noh",value=TRUE),atom.select(res,"noh",value=TRUE))
```

    ##  [1]  0.506  4.310 11.022 10.359  4.781 10.956 10.918  3.704 10.905 10.994
    ## [11] 10.432 10.328 10.846 11.208  8.324  8.935  8.272  8.870

\#\#Normal Mode Analysis of protein structure

``` r
library(bio3d)
pdb <- read.pdb("1HEL")
```

    ##   Note: Accessing on-line PDB file

``` r
modes <- nma(pdb)
```

    ##  Building Hessian...     Done in 0.02 seconds.
    ##  Diagonalizing Hessian...    Done in 0.08 seconds.

``` r
plot(modes, sse=pdb)
```

![](Class13_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Visualize NMA results
mktrj(modes, mode=7, file="nma_7.pdb")
```
