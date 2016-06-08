# CompoundSearch

This repo is meant to help you figure whether a compound of interest is similar to a compound that is already known.
Given a target composition, this code will look through a library of compounds (that you define) and see whether
the target is *similar* to any in that library. The similarity is judged by fuzzy matching techniques of the formula,
which is the specialty of this code. The exact formula does not need to match - we are looking for *similar* formulas.
The greater the similarity, the higher score of the match.

To use this code, run the following as an example:

```python
cs = CompoundSearch()
cs.search("InCl")
```

This code uses the default thermoelectrics library file found in the **compounds_list/thermoelectrics** subfolder.

