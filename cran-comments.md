## R CMD check results

There were no errors, warnings or notes in checking on my local machine (Mac OS).

The package was built on MacOS and also tested using the platform Windows Server 2022, R-release, 32/64 bit (via Rhub) and on winbuilder

I've made a few alterations to the package since the previous manual check with CRAN as follows:

Edits to one examples file to remove a 'commented' example, and ensure examples run in <5 seconds.  Removal of options(warn=-1).
