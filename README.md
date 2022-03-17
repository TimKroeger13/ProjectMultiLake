# ProjectMultiLake
**Name of Plots**

### How to read The Plot names:

1. All Data with CUT_ in there name have there Prior Phase / Pioneer Phase Cut Out. Where the values were cuttet can be seen in the "CorrectionPlots".
2. Data That has Loess in his name is Loess fitted.
3. Data with the name Transform in it is scaled: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/scale
It is both centerd and rescaled.
4. When data is just called mean Data like "Evenness Mean20000" then the line is jsut the mean value in 100 year steps. These means get calculated by using a Ttest.
5. Vector in the name of a dataset indicated that the data is vector based. That means that only the slope between the the interpolated points where used as values and not the values itself. Its therefore relationbased Data.
6. Means are ALWAYS based on the Cutted data!
7. When vector data is transformed the data get first transforment and then the vectors get calculated.
8. When something is called Z_old... its really outdated and its not good using it. Should not be incorrect but not accurate.
