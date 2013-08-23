pyplanet
========

planetary atmosphere code

README file for pyPlanet

*****12/7/17 David DeBoer - created README file
This README file explains the structure and execution of the planetary atmosphere code.  It comprises the following directories:

.(pyModel) : contains the high-level python modules
	planet.py - executive module for automatically running the pipeline
	atmosphere.py - module to compute the atmospheric structure
	microwave.py - module to compute the opacity at each layer
	brightness.py - module to compute the overall brightness temperature
[planetName]: contains planet-specific data input files [Jupiter, Saturn, Uranus, Neptune]
constituents:  sub-directories contain the modules and files to compute the microwave opacity of the constituents
modelSupportInfo:  back-up information/data etc
=================================================================================================================

*****12/12/7 David DeBoer 
Initial posting of pyPlanet (renamed from pyModel) to googleCode
need to set:
	git remote https://code.google.com/p/planetary-atmospheres/
need to add to .netrc
	machine code.google.com/p/planetary-atmospheres/ login username@gmail.com password [currently Vc9aw3ZE5qR6]
then you can:
	git push https://code.google.com/p/planetary-atmospheres/ 
or
	git clone "
Additional files/directories:
	raypath.py - used to calculate raypaths in atmosphere
	regrid.py - regrids data to common array
[planetName]:  added tweakFile.py to contain planet/run dependent atmospheric tweaks
Logs/:  directory where log files are written
Output/: directory where output files are written

*****13/3/20 David DeBoer
git commit -a -m "This skips the adding/staging part..."

