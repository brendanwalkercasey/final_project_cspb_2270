# Wavelet Compression of Imagery
Final Project from CU Boulder's Data Structures 2270, under Applied Computer Science Bachelors Program.


## Overview:
Wavelet compression is a popular method for image compression, which leverages wavelet transforms to represent an image in a more efficient way, typically achieving high compression ratios with minimal loss of quality. The general idea is to decompose the image into a set of wavelet coefficients that can be quantized and encoded more effectively than the original pixel values.  The wavelet transform is a mathematical tool that decomposes an image into different frequency components at multiple scales. These components capture both the low-frequency information (which contains the smooth, overall features of the image) and high-frequency details (which capture sharp edges and fine textures).

In simple terms, it’s akin to applying a filter to the image at different scales, allowing you to see both coarse and fine details in the image. By representing the image in terms of wavelet coefficients, we can selectively retain the most significant information while discarding less important data.

## Details:
- Project data structure implemented:
    - tree (for encoding)
    - bit shifting (for grayscale conversion)
- Short explanation of the data structure:

Using a set of functions, and variables, this process performs the following:

1) Read the TIFF Image: The image is read into a buffer (1D vector) and its dimensions are retrieved.
2) Convert to Numeric Data: The 1D buffer is converted to a 2D array (image) of float values for wavelet transform.
3) Apply Wavelet Transform: The 2D Daubechies wavelet transform is applied in both row and column directions.
4) Post-process and Encode: The transformed data is normalized back to the 0-255 range and saved as a new TIFF image using LZW compression.


Variables:

    Wavelet Decomposition with Daubechies-4 (using low/high pass filters):
        
        The Daubechies-4 filter is used for the 2D wavelet transform. The decomposition into sub-bands happens in two steps: applying the wavelet transform to the rows and then to the columns:

        Low-pass filter (h0, h1, h2, h3)
        high-pass filter (g0, g1, g2, g3)
        
        For each row, the daubechies1D function is called, which applies the Daubechies-4 filter (low-pass and high-pass) to the data. After processing, the data is split into two parts: 1) low frequencies and 2) high frequencies.
        The same process happens for every column in the daubechies2D function, further splitting the image into frequency bands.

        [more needed, tbd]


Functions:

    readTiffImage():
        This function opens a TIFF file, retrieves image dimensions (width and height), as well as the total number of samples per pixel (samples_per_pixel), and reads the pixel data into a buffer (a 1D vector). Values are unsigned 8-bit integers from grayscale images (0-255).  Since my project focuses on implimenting my own image compression, encoding and decoding, this part of the code uses built-in functions from the libTIFF library (https://libtiff.gitlab.io/libtiff/) to open, and read image width, height, and pixel samples (samples per pixel).  
    
- How to run code:
- Write up:


