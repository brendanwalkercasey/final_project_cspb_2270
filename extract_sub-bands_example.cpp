/*
To save each of the 4 sub-bands (LL, LH, HL, HH) after applying the 2D Daubechies-4 wavelet transform, can I modify the code to extract the sub-bands and save them separately as TIFF images?
i.e. Split the Image into Sub-bands
After performing the 2D Daubechies-4 transform (daubechies2D), you need to split the resulting image matrix into the four sub-bands: LL, LH, HL, and HH... and save as separate images?

LL (Low-Low): Contains low frequency components in both the rows and columns
LH (Low-High): Contains low frequency components in the rows, and high frequency in the columns
HL (High-Low): Contains high frequency components in the rows, and low frequency in the columns
HH (High-High): Contains high frequency components in both the rows and columns

