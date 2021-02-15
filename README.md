# holeFiller

Elizabeth B Liddle, Marjie Jansen, Paul M Briley

**Citations:**
  
Briley PM, Liddle EB, Simmonite M, Jansen M, White TP et al. (2020). Regional Brain Correlates of Beta Bursts in Health and Psychosis: A Concurrent Electroencephalography and Functional Magnetic Resonance Imaging Study. Biological Psychiatry: Cognitive Neuroscience and Neuroimaging. Published online ahead of print. doi: 10.1016/j.bpsc.2020.10.018  
  
Liddle EB, Bates AT, Das D, White TP, Groom MJ, Jansen M, et al. (2013). Inefficient cerebral recruitment as a vulnerability marker for schizophrenia. Psychol Med 43:169–182. doi: 10.1017/S0033291712000992  
  
**Installation:**

Ensure the three HoleFiller functions are in your Matlab path (HoleDiagnostics, MultipleImputations_SPMest, CombiningEstimates)

**Diagnostics:**

Call HoleDiagnostics as a function

You can either type HoleDiagnostics into the command window, or you can assign an output variable e.g.:

DodgyFiles=HoleDiagnostics

In which case you will end up with a variable called “DodgyFiles” that will contain the filepaths of the subjects you decide need further investigation, normally those with concerningly small masks as indicated by the bar chart that will be displayed.

You will then be presented with a series of dialog boxes – follow the on-screen instructions.  The first will ask you to navigate to the folder containing your initial second-level analysis.


**Multiple Imputations:**

To impute missing data, again, you will need your original second-level analysis. It is best to have already run all the second-level contrasts you are interested in, as these will then also be output as part of the multiple imputations function.

Then type:

MultipleImputations_SPMest(OriginalDirectory,Threshold,NumberImputations)

where OriginalDirectory is the full path to your original second level analysis; Threshold is either the maximum number or the maximum proportion of missing values you want to impute (if you enter a value >=1, the function will interpret this as a number of voxels (rounding if necessary), and if you enter a proportion of unity, it will interpret this as a proportion of missing values); and NumberImputations is the number of imputations you want to run (5 is standard).

The results of the MI analysis will be saved in a new folder (ImputationsDirectory) that is named as [OriginalDirectory '_imputed']

For each imputation, the function will write a new first level con image (named ImpNo1_con0001.img etc) for the relevant contrast in each subject file, and run a new second level analysis on each subject, which will be stored in the ImputationsDirectory in a separate subdirectory (ImpNo1, ImpNo2 etc). 

When these have been done, the results will be combined using CombiningEstimates.m, and written to an additional subdirectory, called MeanAnalysis which will become the new working directory. You can then call SPM and interrogate this analysis, which should contain all the contrasts specified in your original second level analysis.
