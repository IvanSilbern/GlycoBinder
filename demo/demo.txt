Demonstration data set

As a test data set, we provide an IgM_TMT0.raw file. It is a tryptic digest of a purified IgM sample labeled with TMT0 reagent.
The file is located in the “demo” folder together with a Human_IgM.FASTA file containing amino acid sequences of the two human proteins, IgM and IgJ, respecitvely.
To test the performance of the GlycoBinder, download the contents of the “demo” folder (e.g. into C:/data/Glycobinder/demo),
copy the current version of GlycoBinder into it and execute in the command line using following parameters:

C:/data/Glycobinder/demo>Rscript.exe "GlycoBinder.R" --wd "C:/data/Glycobinder/demo" --reporter_ion TMT0 --no_second_search

If you download files from GitHub using git bash, please first install git lfs (https://git-lfs.github.com/) that is aimed at handling large files (e.g. the example raw file). If you use web interface for downloading (“download ZIP”), you will download a placeholder for IgM_TMT0.raw file. To download the actual file, find it in the GitHub repository and click on “View raw”. Save the file in your local “demo” folder.

The execution takes around 5 min on a desktop computer running Windows 10 and equipped with Intel Core i7-6700 CPU (64 bit) and 32 Gb of RAM.
Beware that the execution time will scale up with the complexity of the data set provided.

GlycoBinder output is located within pglyco_output folder. There are 67 glycoforms (pGlyco_glycoforms.txt), 45 unique glycan compositions (pGlyco_glycans.txt) and 4 glycosylation sites (pGlyco_glycosites.txt) identified in the data.
