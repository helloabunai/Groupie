Groupie: Batching sequence repeat count distributions
=========================================================
Simple script for colleagues to use for batching repeat counts into specified groups.
Useful in low read-count situations, for further data analysis.


Installation Prerequisites
==========================

It is best practice to ensure your version of 'setuptools' is up-to-date before trying anything. If you are not using
an administrator account currently, you may need to use sudo.

    $ easy_install -U setuptools

And for proper package installation:

    $ cd ~/src/
    $ python setup.py install

Usage
=====

Here's how to use it for now. To launch the program when installed, launch via terminal:

    $ groupie -flags

All flags for the application are:

    -sam :: A path to a sam file, or a folder of sam files as input
	-rng :: A numerical range to batch reads into [Default is 0-10]
    -cpu :: Number of CPU threads to utilise [Default is maximum]
    -out :: A path to a desired output folder [Default is $HOME/Groupie/out/]
