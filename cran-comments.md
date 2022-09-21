## R CMD check results

There were no errors, warnings or notes in checking on my local machine (Mac OS).

The package was built on MacOS and also tested using the platform Windows Server 2022, R-release, 32/64 bit (via Rhub) and on winbuilder

I've made a few alterations to the package since the previous manual check with CRAN as follows:

Edits to the title and description fields of the description, edits to functions so they would not print to console unnecessarily, edits to help files (so that the return value is listed for all functions), small changes to examples to ensure all examples run, are not commented out, and the replacement of the dontrun directive with donttest in the examples.
