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

1) Read the TIFF image into a buffer: The image is read into a buffer (1D vector) and its dimensions are retrieved.
2) Apply wavelet compression (including the Daubechies transform and quantization).
    2 a) Convert to Numeric Data: The 1D buffer is converted to a 2D array (image) of float values for wavelet transform.
    2 b) Apply Wavelet Transform: The 2D Daubechies wavelet transform is applied in both row and column directions.
3) Write the compressed image to a new TIFF file.
4) Post-process and Encode: 
    4 a) The transformed data is normalized back to the 0-255 range and saved as a new TIFF image using LZW compression.
    4 b) Reconstruct the image using the inverse wavelet transform.
    4 c) Write the reconstructed image to a new file.




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
    
## How to run code:
## Write up:

Description:

This code performes the following:

1) Wavelet Decomposition (Daubechies-4):
    - A Daubechies-4 (Db4) wavelet transform, which decomposes the image into frequency components (low and high frequencies) at multiple scales.
    - The forward transform is performed in two steps:
        - Row-wise: The 1D Daubechies-4 transform is applied to each row.
        - Column-wise: The same transform is applied to each column of the image.
    - This step decomposes the image into four sub-bands:   
        - LL (Low-Low): Low frequency in both rows and columns.
        - LH (Low-High): Low frequency in rows, high frequency in columns (horizontal detail).
        - HL (High-Low): High frequency in rows, low frequency in columns (vertical detail).
        - HH (High-High): High frequency in both rows and columns (diagonal detail).

2) Quantization:
This project makes use of Bit-Shifting Quantization: 
    - In the applyBitShiftingQuantization() function, quantization is applied by shifting the wavelet coefficients. The goal is reduced precision for each applicable coefficient (through right-shifting), which helps in reduce the size of the data (thus, achieving compression). This technique is a form of lossy compression, which discards less significant information (small coefficients) and keeps the more important ones (large coefficents).


3) Compression:
    - Compression: After the wavelet transform is applied, quantization is performed on the transformed coefficients, and then normalization of the data back to the 0-255 range (suitable for image storage).
    - Reconstruction: inverseWaveletTransform() function inversely transforms the image, reconstructing it from the quantized coefficients. This step involves reversing the wavelet transform.  After applying the inverse transform, the image is reconstructed and written to a new TIFF file.  This technique results in some loss of image quality due to the quantization.


## Analysis

Bit Shift Values:
The bitShiftAmount is used for applying bit-shifting quantization, which reduces the precision of image data by right-shifting the pixel values. This shift effectively narrows the range of pixel values, making it useful for compression by reducing the amount of data required to represent the image. As expected, higher bitShiftAmount values (e.g., 4 or 5) result in greater compression, but at the cost of reduced image quality due to more aggressive quantization. This can lead to a noticeable loss of detail, particularly in areas with fine textures or smooth gradients. On the other hand, lower bitShiftAmount values (e.g., 1 or 2) preserve more detail but achieve less compression. While this may not significantly reduce data storage, it still offers some level of quantization to decrease file size.  This can be demonstrated below:

    - bitShiftAmount = 1, or 2 (Light Compression: preserve image quality while still reducing the file size, useful for noise reduction)
    ![Example Image](v11_shift2_output_compressed.tif)
    - bitShiftAmount = 3 (Moderate Compression: balances file size and image quality)
    - bitShiftAmount = 4, or 5 (Heavy Compression; reduce file size, more important than preserving image quality)

