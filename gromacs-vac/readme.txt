#################################################################
Title     : GROMACS vacuum simulation scripts for version 4.6.7
Author    : ycshi @ CSRC
TimeStamp : 2019.6
#################################################################

I. Start a project
    Make a folder, for example "~/Documents/startup", which should contain input files such as "xxx.mdp", "xxx.gro" and "xxx.top". Then go to a place where you want to create this project, run
    '''
    $ ./initiate.sh ~/Documents/startup
    '''

II. Run vacuum simulation
    The simplest way to start simulation is to run
    '''
    $ ./timer.sh
    '''
    The scripts will generate input files, submit jobs and monitor job status automatically.
    Each simulation lasts 1 nanosecond, then output ".gro" & ".top" file will be checked and updated to generate new input files for next 1ns simulation.
    If you just want to submit only next 1ns simulation, run
    '''
    $ ./autoload.sh
    '''
    Other bash and python scripts are supported files.

III. Analysis
    There is a "analysis" folder in the project for further analyzing.
    All of the analysis scripts in "mdana" project is encouraged to be used in your program.
