26 October 2011

GPT: General Particle Tracer version 3.10
Pulsar Physics
info@pulsar.nl
www.pulsar.nl/gpt

This README file contains errata and information not covered in
the GPT Manuals, as well as additional troubleshooting tips.

KNOWN INSTALLATION ERRORS:
- Impossible to install side-by-side with previous GPT releases

KNOWN ERRORS:
- GDF2DXF produces 2D output incompatible with AutoCAD
- Font is not stored in plot templates

ADDITIONAL INFORMATION:
- The latest news about GPT can be found at: 
  http://www.pulsar.nl/gpt/news.html

- 32-bit Custom Elements rely on the Microsoft Visual C++ 2008 'Express' Edition.
  It can be downloaded free of charge from:
  http://www.microsoft.com/express/download/

  64-bit Custom Elements require the 'Standard' or 'Education' version.
  
  To add custom elements with Microsoft Vista, please run GPTwin with
  Administrative permissions (Run As Administrator).

NEW IN GPT VERSION 3.10
- New elements: TM110cylcavity, map3D_TM, map3D_Ecomplex, map3D_Hcomplex
    setellipse, setscale, scatterbitmap, pointchargeset
- Multicore spacecharge3Dtree
- Compiler switch to Microsoft Visual C++ 2010
- Perfect screens without interpolation errors
- Snapshots

NEW IN GPT VERSION 3.03
- Mixed charge capabilities for spacecharge3Dmesh

NEW IN GPT VERSION 3.02
- 'No coarser grid error' fix in spacecharge3Dmesh
- Spacecharge3Dmesh solves in rotated frame
- Multi-threaded MR

NEW IN GPT VERSION 3.01
- Compiler switch to Microsoft Visual C++ 2008
- Native 64-bit versions of all MS-Windows executables

NEW IN GPT VERSION 3.00
- spacecharge3Dtree

NEW IN GPT VERSION 2.82
- Performance increase starting >1M particles
- Native 64-bit versions of GPT and GDFA
- Export (density)data
- If statements in inputfile
- HTML-Help for all GPTwin commands

NEW IN GPT VERSION 2.81
- Export data commands
- Corrected setfile for nested groups
- Individual phase-space projections can be read from file
- Rewritten screen element optimized for continuous beams
- Option for external color-table in color-density plots
- Corrected magdipole
- Corrected job control under Vista

NEW IN GPT VERSION 2.80
- New windows installer with 'repair' functionality
- GDF2GDF program to combine GDF files
- Automatic timestep reduction to prevents skipping small local elements
- New elements: rectmagnet, sectormagnet, setxydistbmp, ehole
- All keywords and variables are parsed case-insensitive
- Several minor bugfixes

NEW IN GPT VERSION 2.71
- Colored scatterplots
- Requires free Microsoft Visual Toolkit 2003 compiler
- Improved error detection and recovery in GDFA
- Various small changes in MR such as negative stepsizes

NEW IN GPT VERSION 2.70
- 3D mesh-based space-charge model (spacecharge3Dmesh)
- Animation support in GPTwin
- Logarithmic scaling possible
- Nested MR with single inputfile
- MPI version of MR
- Export to .bmp and .avi
- FEL simulations (gauss00mf)
- New elements: setxyzgrid, gauss00mf, setcurvature, setstartxyzgrid, collision, setrmacrodist,
    TErectcavity, TMrectcavity, Tm010cylcavity, magdipole, stdxyzmax, Gminmax
- GDFdemo program added
- 1D GDF based field-maps

NEW IN GPT VERSION 2.61
- GDF2SDDS
- New elements: scatter

NEW IN GPT VERSION 2.60
- Stop button
- Huge file (>4GB) support
- Start particles as function of time with settdist
- Hammersley initial particle distributions
- Descriptive window titles
- New elements: settdist
