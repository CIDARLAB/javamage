# Javamage

## What is Javamage

Javamage is the Java part of Merlin

Merlin consists of ...

* Front-End: This is how users can interact with Merlin inputs and outputs 
* Server : This handles simple GET and POST requests made by the Front End.  The server then uses Javamage to get run optimization
* Javamage : This is a Java package for performing Merlin Analysis on the package

## How to use Javamage

See the javadoc for a more detailed explanation of how everything works.

###  Setup and Requirements

Merlin has the following Java based requirements. These are .jar files that can be found in the `/javamage/lib` folder

* BioJava 3.0 [Sequenced, Alignment, Core]
* Javaplot [A plotting library for anaylsis]

Merlin has requires the follow EXTERNAL dependencies

* `makeblastdb`:  , comes with NCBI Blast+ Standalone
* `blastn`:  , comes with NCBI Blast+ Standalone
* `MFOLD` : This requires UNAFOLD, available at `http://mfold.rna.albany.edu/?q=DINAMelt/software`

after installing these three dependencies.

To configure javamage to find these binaries, write down in the config.txt file found in the working directory the process name followed by the path to the binary

You can use the `which` command to find out where these binaries live.

Below is a sample config.txt entry

<pre>
makeblastdb /usr/local/ncbi/blast/bin/makeblastdb
blastn /usr/local/ncbi/blast/bin/blastn
MFOLD /usr/local/bin/hybrid-ss-min
</pre>

> ensure that you follow the order of this format

###  Getting Started

Look at the `javamage/test/Unit/Merlin.java` file to get an example of how to use java ma

Open the file and change the 'merlin test dir' string  to link to YOUR wordking directory folder.

After this setup the correct file names. These files should be in the working directory folder

* `genomeFileName` : A Fast file of the genome
* `targetFileName` : An optMAGE target File name
* `parametersFileNamie` : An optMAGE parameter File

Using these 4 parameters we create a new Merlin object

We enable and disable certain procedures like plotting.

We can set Switches, like the BLAST Scoring method, with the `mage.Switches.Blast.method=2` line.

The finally the `.optimize()` method will run the Merlin Optimization analysis.
Then on the last commented line is allows you to compare merlin results to optmage results that are provided in an OUTOligos file in the resources folder)

```
Java Analysis Tool with optMAGEv0.9
Created by Samir Ahmed 2011-2012
```
