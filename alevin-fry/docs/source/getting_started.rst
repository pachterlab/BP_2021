Getting Started
===============

There are a (growing) number of different sub-commands: generate-permit-list, collate, and quant. Each of these is invoked as a command passed as the first argument to the alevin-fry executable. For example, to run the generate-permit-list command, one would run:

.. code:: bash

    $ alevin-fry generate-permit-list --help

This should then show the following:

.. code:: bash

    $ alevin-fry generate-permit-list --help
    alevin-fry-generate-permit-list 0.0.1
    Avi Srivastava, Rob Patro
    Generate a permit list of barcodes from a RAD file

    USAGE:
        alevin-fry generate-permit-list [FLAGS] --input <input> --output-dir <output-dir> --expect-cells <expect-cells> --force-cells <force-cells> --valid-bc <valid-bc>

    FLAGS:
        -h, --help             Prints help information
        -k, --knee-distance    attempt to determine the number of barcodes to keep using the knee distance method
        -V, --version          Prints version information

    OPTIONS:
        -e, --expect-cells <expect-cells>    defines the expected number of cells to use in determining the (read, not UMI) based cutoff
        -f, --force-cells <force-cells>      select the top-k most-frequent barcodes, based on read count, as valid (true)
        -i, --input <input>                  input RAD file
        -o, --output-dir <output-dir>        output directory
        -b, --valid-bc <valid-bc>            uses true barcode collected from a provided file


Running the alevin-fry pipeline
-------------------------------

First, we need to generate the RAD file using alevin.  For a chromium v2 set of read files, the command would look like the following:

.. code:: bash
    $ salmon alevin -lISR --chromium -1 <read1_files> -2 <read2_files> -o <alevin_odir> -i <index> -p <num_threads> --tgMap <tg_map> --justAlign --sketchMode 

Given the output directory generated by alevin above, the next step is to let alevin-fry generate the permit list.  Here we use the "knee" method `-k`.

.. code:: bash 
    $ alevin-fry generate-permit-list --input <alevin_odir> --expected-ori fw --output-dir <fry_odir> -k

Next, given the permit list and barcode mapping (which resides in the `<fry_odir>` directory), we collate the original RAD file using the command below.

.. code:: bash 
    $ alevin-fry collate -i <fry_odir> -r <alevin_odir> -t <num_threads>

Finally, we quantify the collated rad file using the `cr-like` resolution strategy using the `quant` command below.

.. code:: bash 
    $ alevin-fry quant -i <fry_odir> -m <tg_map> -t <num_threads> -r cr-like -o <fry_odir> 

Note that with the exception of the `generate-permit-list` command, the other `alevin-fry` commands are designed to scale well with the number of provided threads. Thus, if you have multiple threads to use, then you can provide the appropriate argument to the `-t` option.