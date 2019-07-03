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

II. Edit submit.sh
    Here "submit.sh" should be re-edit to submit jobs in your environment.
    The default script is designed for Sugon cluster.

III. Run vacuum simulation
    The simplest way to start simulation is running
    '''
    $ ./timer.sh
    '''
    The scripts will automaticaly generate input files, submit jobs and monitor status.
    Each simulation lasts for 1 nanosecond. The output ".gro" & ".top" files will be updated to generate new input files for the next 1ns' simulation.
    If you just want to submit the only next 1ns simulation, run
    '''
    $ ./autoload.sh
    '''
    Other bash and python scripts are for supporting.

III. Analysis
    There is a "analysis" folder in the project for further analyzing.
    All of the analysis scripts in "mdana" project is encouraged to be used in your program.
