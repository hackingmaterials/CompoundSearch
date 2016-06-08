# CompoundSearch

This repo is meant to help you figure whether a compound of interest is similar to a compound that is already known. Given a target composition, this code will look through a library of compounds (that you define) and see whether the target is *similar* to any in that library. The similarity is judged by fuzzy matching techniques, which is the (intended) specialty of this code - although keep in mind that these are just estimates. The exact formula does *not* need to match - we are looking for *similar* formulas. The greater the similarity, the higher **score** of the match. The convention is that an exact match to the formula has a **score** of 100.

To use this code, run the following as an example:

```python
cs = CompoundSearch()
cs.search("Pb2TeSe")
```
Note that if you run "compoundsearch.py", the main method at the end will essentially do this.

which will return something like this (first couple of output lines):
```
formula_pretty	score	matches	line_no	tag1	tag2	tag3
FeCuO2	15	['anonymized formula']	5	 Gonçalves, A. P. & Godart, C. New promising bulk thermoelectrics: Intermetallics, pnictides and chalcogenides. Eur. Phys. J. B 87, 42 (2014).	p2	
BiSb	15	['fuzzy formula']	22	 Gonçalves, A. P. & Godart, C. New promising bulk thermoelectrics: Intermetallics, pnictides and chalcogenides. Eur. Phys. J. B 87, 42 (2014).	p7	
AgSbTe2	15	['anonymized formula']	40	 Gonçalves, A. P. & Godart, C. New promising bulk thermoelectrics: Intermetallics, pnictides and chalcogenides. Eur. Phys. J. B 87, 42 (2014).	p9	
GaAgTe2	15	['anonymized formula']	41	 Gonçalves, A. P. & Godart, C. New promising bulk thermoelectrics: Intermetallics, pnictides and chalcogenides. Eur. Phys. J. B 87, 42 (2014).	p9	
GaCuTe2	15	['anonymized formula']	42	 Gonçalves, A. P. & Godart, C. New promising bulk thermoelectrics: Intermetallics, pnictides and chalcogenides. Eur. Phys. J. B 87, 42 (2014).	p9	
AgSbTe2	15	['anonymized formula']	47	 Gonçalves, A. P. & Godart, C. New promising bulk thermoelectrics: Intermetallics, pnictides and chalcogenides. Eur. Phys. J. B 87, 42 (2014).	p10	
GeTe	30	['group formula']	48	 Gonçalves, A. P. & Godart, C. New promising bulk thermoelectrics: Intermetallics, pnictides and chalcogenides. Eur. Phys. J. B 87, 42 (2014).	p10	
SnTe	45	['group formula', 'fuzzy formula']	51	 Gonçalves, A. P. & Godart, C. New promising bulk thermoelectrics: Intermetallics, pnictides and chalcogenides. Eur. Phys. J. B 87, 42 (2014).	p10	
TePb	45	['group formula', 'fuzzy formula']	52	 Gonçalves, A. P. & Godart, C. New promising bulk thermoelectrics: Intermetallics, pnictides and chalcogenides. Eur. Phys. J. B 87, 42 (2014).	p10	
```

The highest score of these is for ``SnTe`` and ``TePb`` since they have a score of 45 (maximum amongst the rows listed). The ``GeTe`` is close in terms of a match to the target.

This code uses the default thermoelectrics library file found in *compounds_list/thermoelectrics.txt*. This library file contains some known thermoelectrics compositions; the code is searching to see which of the known compositions best matches the query "InCl". There are no exact matches but a few formulas that have at least passing similarity in some sense.

You can modify or create a new library file as desired. The library file format is simple. Each line contains a formula that is searched. If there is a line that begins with "#", "##", or "###", the line used to tag the entries below with a 1st level, 2nd-level, or 3rd-level tag (respectively). The tags will be reported along with the compound match and can be used to help track down where you found the compound. For example, you can use the 1st level tag to give a paper reference, a 2nd level tag to give the page number, and a 3rd level tag for a note. This is the convention used in the default library file.

## FAQS

### What algorithm(s) are being used to judge formula similarity?

There are a few different types of matches we are looking for:
* exact formula matches
* same chemical system
* same anonymized formula
* formula match if elements are replaced by their group number
* formula match if elements by placeholders such as "A" for cation, "Tm" for transition metal, "X" for anion, etc.

Each type of match has a different score that is added if the match is confirmed. Suggestions for algorithm improvements are welcome.

### I am getting a lot of results reported back that don't look very similar.
This is common - you basically need to decide for yourself which formulas are truly similar and which ones are not from the list provided. The code is only there to feed you ideas for what known compounds might be similar to your composition. It is not a complete solution by any means.

### Can you improve the interface? e.g., allow me to sort the results by score? Or give me a command line interface?

If there is enough interest in this code, these things are possible.


