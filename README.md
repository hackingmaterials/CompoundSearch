# CompoundSearch

This repo is meant to help you figure whether a compound of interest is similar to a compound that is already known. Given a target composition, this code will look through a library of compounds (that you define) and see whether the target is *similar* to any in that library. The similarity is judged by fuzzy matching techniques, which is the specialty of this code. The exact formula does not need to match - we are looking for *similar* formulas. The greater the similarity, the higher **score** of the match. The convention is that an exact match to the formula has a **score** of 100.

To use this code, run the following as an example:

```python
cs = CompoundSearch()
cs.search("InCl")
```

This code uses the default thermoelectrics library file found in *compounds_list/thermoelectrics.txt*. This library file contains some known thermoelectrics compositions; the code is searching to see which of the known compositions best matches the query "InCl". There are no exact matches but a few formulas that have at least passing similarity in some sense.

You can modify or create a new library file as desired. The library file format is simple. Each line contains a formula that is searched. If there is a line that begins with "#", "##", or "###", the line used to tag the entries below with a 1st level, 2nd-level, or 3rd-level tag (respectively). The tags will be reported along with the compound match and can be used to help track down where you found the compound. For example, you can use the 1st level tag to give a paper reference, a 2nd level tag to give the page number, and a 3rd level tag for a note. This is the convention used in the default library file.

