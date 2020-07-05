# autotc
This is a (not particularly well thought out) set of scripts that I use to generate lots of TeraChem jobs with slightly different parameters, and/or for many different structures. There is also a very barebones process scheduler included, though the main script can use any process scheduler (or none at all) if the user prefers. 



Functionality is as follows:
There are two main configuration files, job\_params.ini and xyz\_params.ini, and a job will be created for each combination of sections in job\_params.ini and xyz\_params.ini. Entries in job\_params.ini have the format "[terachem keyword]=[alias]", and entries in xyz\_params.ini have the format "[alias]=[argument to terachem keyword]". This is probably best explained with an example:


```ini
xyz_params.ini:
----

[ethane]
coords=ethane.xyz
c=4
a=4
[methane]
coords=methane.xyz
c=3
a=4

----
job_params.ini
----

[fomo_casci]
template=fomo_casci.in
closed=c
active=a
[hf_casci]
template=hf_casci.in
closed=c
active=a

```

Note that coords and template are special arguments that refer to files in the templates directory. This will create four jobs, a (4,4)fomo-casci for ethane and methane, and a (4,4)hf-casci for ethane and methane. The directories will be created, and the input templates will be search-and-replaced with the specified terachem keywords and arguments to create the input files. This seems like a convoluted way to set up four jobs, but this can save a lot of time when one needs to benchmark a method on a large set of molecules that all need different active space parameters. In addition, it only takes a couple lines of python to generate a pair of ini files to cover combinations of several varying parameters. Example scripts are in the cfg directory.

