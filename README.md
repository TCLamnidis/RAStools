# RAStools
A toolset to calculate Rare Allele Sharing (RAS).

The input files should be in FreqSum format. The FreqSum format and various tools for working with this format are documented in https://github.com/stschiff/rarecoal-tools.

# Using RASCalculator
* The input FreqSum file can be input using the option `-I`, or as stdin by 
	directly piping the FreqSum to `RASCalculator.py` and omitting this 
	option.
* The output file is defined using the `-O` option.
* The outgroup with which all alleles should be polarised is defined with the 
	`--outgroup` option. This is also the outgroup used for F3 
	calculations. We suggest using Chimp in most cases.
* The maximum allele frequency for the non-Reference allele in can be changed 
	with the `-M` option (default is 10).
* By default, the minimum allele frequency of the non-reference allele is 2, 
	but it can be changed with option `-m` to analyse other slices of RAS.
* The `--ascertainIn` option can be used to define the comma-separated set of 
	populations that should be used for rare allele ascertainment. 
	By default, all non-Left populations are used to compute allele counts. 
* The `--calculateFor` option can be used to define the comma-separated set of 
	populations that should be tested for RAS with each Left population. 
	By default, this will include all non-Left populations.
* A set of "Left" populations can be defined as comma-separated list and 
	provided to the `--LeftPops` argument. 
* A cutoff for missingness in the ascertainment populations can be added by 
	using the `--MissingnessCutoff` option (default=0.0). A cutoff of `0.1` 
	will exclude variants where more than 10% of the individuals in the 
	ascertainment populations contain a missing call.
* With the `-P` flag it is possible to calculate RAS on shared variants that are
	privately shared between a Left population and a tested population.
* The `-NT` flag will exclude variants that are transitions from the RAS 
	calculation. Especially useful for calulating RAS with ancient 
	individuals, where DNA damage may affect the results.
* For use with non-human data, the `--NrChroms` option can be used to change the
	number of chromosomes.
* By default, `RASCalculator.py` will print out RAS for all variants between the
	specified allele counts. 
	The `--details` option will enable printing of additional RAS 
	calculations for each separate allele count.
* Finally, the `--f3FreqCutoff` option can be used to change the cutoff for 
	minor allele frequency in the ascertainment populations, used for F3 
	calculation. The default value is 0.05, which we have found to reduce 
	sensitivity on ancient DNA damage, and exclude rare variants from the 
	calculation.
* The above information can also be found in the command line with the option `-h`.
```
usage: RASCalculator.py [-h] [-I <INPUT FILE>] [-O <OUTPUT FILE>] -o POP
                        [-M <MAX ALLELE COUNT>] [-m <MIN ALLELE COUNT>]
                        [-a POP1,POP2,...] -c POP1,POP2,... -L POP1,POP2,...
                        [-x <CUTOFF>] [-NT] [-P] [-C <INT>] [-d]
                        [--f3FreqCutoff <FREQ>]

Compute rare allele sharing statistics between two populations with respect to
an outgroup, as well as outgroup F3 statistics. Also preforms error estimation
using jackknifing, using the number of observed sites for normalisation.

optional arguments:
  -h, --help            show this help message and exit
  -I <INPUT FILE>, --Input <INPUT FILE>
                        The input freqsum file. Omit to read from stdin.
  -O <OUTPUT FILE>, --Output <OUTPUT FILE>
                        The output file. Omit to print in stdout.
  -o POP, --outgroup POP
                        The outgroup to polarize all alleles with. As a useful
                        standard choice, use Chimp as outgroup
  -M <MAX ALLELE COUNT>, --maxAF <MAX ALLELE COUNT>
                        The maximum number of alleles (total) in the reference
                        populations. The default maximum allele value is 10.
  -m <MIN ALLELE COUNT>, --minAF <MIN ALLELE COUNT>
                        The minimum number of alleles (total) in the reference
                        populations. The default minimum allele count is 2.
  -a POP1,POP2,..., --ascertainIn POP1,POP2,...
                        The populations to ascertain in. Rare alleles will be
                        defined as alleles with minor allele freuency between
                        <minAF> and <maxAF> in this set of populations. By
                        default, all non-Left populations will be used for
                        ascertainment.
  -c POP1,POP2,..., --calculateFor POP1,POP2,...
                        The populations to be used in RAS calculations. RAS
                        will be calculated between each Left population and
                        each of these populations. By default, RAS will be
                        calculated for all non-Left populations.
  -L POP1,POP2,..., --LeftPops POP1,POP2,...
                        Set the Left populations/individuals. RAS will be 
			calculated between each Left and all populations in 
			the calculation list.
  -x <CUTOFF>, --MissingnessCutoff <CUTOFF>
                        Missingness cutoff proportion for Right populations.
                        E.g. 0.1: If more than 10% of individuals in Right
                        populations show missing data, the variant will be
                        ignored. [default=0]
  -NT, --NoTransitions  When present, No Transitions are included in the
                        output. Useful for ancient samples with damaged DNA.
  -P, --Private         Restrict the RAS calculation to privately shared rare
                        variants only.
  -C <INT>, --NrChroms <INT>
                        The number of chromosomes in the dataset. [22]
  -d, --details         Print RAS calculations for each allele frequency, in
                        addition to the total.
  --f3FreqCutoff <FREQ>
                        minimum minor allele frequency for f3 values. Defaults
                        to 0.05. This is recommended to skip rare mutations
                        and reduce the sensitivity of the F3 statistics on
                        sequencing errors and DNA damage.
```
# RAS Output Format

The first two lines of a RAS output start with a `#` and contain the populations and sizes of the populations in the input FreqSum, and the population information about that run. 

    #Available populations in Input File and their respective sizes: {'PopA': 8, 'PopB': 14, 'PopC': 14, 'PopD': 28, 'PopE': 20, 'PopF': 24, 'PopG': 24, 'Ind1': 2, 'Ind2': 2, 'Ind3': 2, 'Ind4': 2, 'Ind5': 2, 'Chimp': 2}
    #Left Populations:  Ind1 Ind2 Ind3
    #Test Populations:  PopA PopB PopC PopD
    #Populations considered for allele frequency calculation (Rights):	PopA	PopB	PopD	PopE	PopF
    #Outgroup: Chimp

This information is followed by the result tables.
The first line of the result table contains the labels for each result column. 
After this line the script will output a table for each population, containing the above information across all allele frequencies as well as one the outgroup F3 statistic.

    TestPop	LeftPop		RAS	    RAS/Mb	              Jackknife Estimator     Jackknife Error   Allele Frequency
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Total [2,10]
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Outgroup F3
    PopB    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Total [2,10]
    PopB    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Outgroup F3
    PopC    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Total [2,10]
    PopC    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Outgroup F3
    PopD    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Total [2,10]
    PopD    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Outgroup F3
    PopA    Individual2 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Total [2,10]
    PopA    Individual2 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Outgroup F3
    PopB    Individual2 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Total [2,10]
    PopB    Individual2 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Outgroup F3
                            ...

If the `--details` option is provided the output will also contain the calculations per allele frequency.

    TestPop	LeftPop		RAS	    RAS/Mb	              Jackknife Estimator     Jackknife Error   Allele Frequency
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	2
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	3
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	4
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	5
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	6
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	7
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	8
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	9
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	10
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Total [2,10]
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Outgroup F3
    PopB    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	2
    PopB    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	3
                        ...

