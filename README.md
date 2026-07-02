# mixed_ll and mixed_nul

## External dependencies
In order to properly install and run `mixed_ll` the following packages must be available:

+ [LHAPDF][2]
+ [OpenLoops][3]

## OpenLoops installation and configuration
You can download `OpenLoops` from git via

    git clone https://gitlab.com/openloops/OpenLoops.git

or, if you are interested in the latest beta version, via

    git clone -b public_beta https://gitlab.com/openloops/OpenLoops.git

and then compile it using the following procedure:

    cd OpenLoops
    ./scons

In order to be able to access all the process libraries needed in this project you need to create an 'openloop.cfg' file adding the following:

    [OpenLoops]
    process_repositories = public, public_beta, nested_soft_coll

You can download and install all the relevant process libraries going into your `OpenLoops` folder and then typing

    ./openloops libinstall pplla
    ./openloops libinstall ppllj_ew
    ./openloops libinstall pplljj_qcdew

[2]: https://lhapdf.hepforge.org/ "LHAPDF"
[3]: https://openloops.hepforge.org/ "OpenLoops"
