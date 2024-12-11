# Wavelet Compression of Imagery
Final Project from CU Boulder's Data Structures 2270, under Applied Computer Science Bachelors Program.


## Overview:
Wavelet compression is a popular method for image compression, which leverages wavelet transforms to represent an image in a more efficient way, typically achieving high compression ratios with minimal loss of quality. The general idea is to decompose the image into a set of wavelet coefficients that can be quantized and encoded more effectively than the original pixel values.  The wavelet transform is a mathematical tool that decomposes an image into different frequency components at multiple scales. These components capture both the low-frequency information (which contains the smooth, overall features of the image) and high-frequency details (which capture sharp edges and fine textures).

In simple terms, it’s akin to applying a filter to the image at different scales, allowing you to see both coarse and fine details in the image. By representing the image in terms of wavelet coefficients, we can selectively retain the most significant information while discarding less important data.

This project uses a Daubechies Wavlet (Daubechies-4): for decomposition of an image into frequency components, bit-shifting quantization: to quantize large values, while smaller coefficients are thresholded out with a sigma value, an inverse wavelet transform: to reconstruct the image from its compressed components, and built-in methods (from libTIFF library) to encode the image after quantization (with Limpel-Ziv-Welch, or LZW Compression).  


## Details:
Data structures and techniques implemented:
    - Vectors (main base structure for storing image data, in the form of 1-D and 2-D arrays, typically with 8-bit, 16-bit, unsigned integers, or float types)
    - Buffers (temporary structures used during 2D Daubechies tranform to store intermediate results)
    - Pointers (namely pointers to TIFF objects, to interact with the libTIFF library in a memory efficent manner)
    - Standard Library Algorithms (standard built-in methods to better organize code, and save time e.g. min_element and max_element for finding min/max values)
    - Bit-Shifting (simple and fast operations with minimal computations to reduce data on individual pixel values, effectively compressing image file with fewer bits)


Using a set of functions, and variables, this code performs the following:
    1) Read the TIFF Image into a Buffer:
        Action: The image is read into a buffer (vector<u8>), where each pixel is stored as a grayscale value (an 8-bit unsigned integer).
        Details: The image's dimensions (width, height) and the number of samples per pixel (e.g., grayscale = 1) are retrieved from the TIFF metadata.
    2) Apply Wavelet Compression (Daubechies-4 Transform):
        2a) Convert to Numeric Data:
            - The 1D buffer (which holds pixel values) is converted to a 2D matrix (vector<vector<float>>) where each value is a floating-point representation of the pixel. This conversion is necessary for performing the wavelet transformation.
        2b) Apply Wavelet Transform:
            - The 2D Daubechies-4 wavelet transform is applied to the image in two steps:
                1) 1D transform is applied to each row of the image.
                2) 1D transform is applied to each column of the image.
    3) Quantization:
        Action: The wavelet coefficients are quantized using bit-shifting. This process reduces the precision of the coefficients for compression.
        Details: Small coefficients are thresholded to zero based on a sigma value (applyBitShiftingQuantization()). Larger coefficients are shifted by a specified amount (e.g., by 2 bits), which reduces their range and thus compresses the data.
    4) Reconstruction (Inverse Wavelet Transform):
        Action: The inverse Daubechies-4 wavelet transform is applied to reconstruct the image from its wavelet sub-bands.
        Details: The image is reconstructed by applying the inverse wavelet transform (first on columns, then on rows), reversing the effect of the initial wavelet transform.
    5) Post-process and Encode:
        5a) Normalize the Data: The transformed image data is normalized back to the 0-255 grayscale range for proper display. This normalization ensures that the pixel values fit within the standard 8-bit range.
        5b) Save the Processed Image: The compressed (and potentially quantized) image is saved to a new TIFF file using LZW compression for file size reduction.5c) Reconstruct and Save the Image: The image is reconstructed using the inverse wavelet transform (as described earlier) and then saved to another TIFF file as the "reconstructed" version.
    6) Additional Steps: 
        - Sub-Band Decomposition: After applying the wavelet transform, four sub-bands are created to illustrate details of directional components of the image:
            - LL (approximation) (visible in top left of compressed output tiff)
            - LH (horizontal details) (visible in bottom left of compressed output tiff)
            - HL (vertical details) (visible in top right of compressed output tiff)
            - HH (diagonal details) (visible in bottom right of compressed output tiff)
        While these details are visible in the compressed file (see list above), the code for generating individual output tiff files (e.g. LL tiff, or HH tiff) did not correctly split and resize the original compressed output image.  Therefore, generating individual tiff files for each sub-band was left as optional under the main() function.  

Variables:

    Wavelet Decomposition with Daubechies-4 (using low/high pass filters):
    
        The Daubechies-4 filter is used for the 2D wavelet transform. The decomposition into sub-bands happens in two steps: applying the wavelet transform to the rows and then to the columns:

        Low-pass filter (h0, h1, h2, h3)
        high-pass filter (g0, g1, g2, g3)
        
        For each row, the daubechies1D function is called, which applies the Daubechies-4 filter (low-pass and high-pass) to the data. After processing, the data is split into two parts: 1) low frequencies and 2) high frequencies.
        The same process happens for every column in the daubechies2D function, further splitting the image into frequency bands.

    Saving Daubechies-4 to separate files: 

        Initialize minVal and maxVal: The main reason for using constants is for defining initial extremes for the minVal and maxVal variables. When starting with largest possible, and smallest possible float values, we ensure that every value in the collection will either be smaller than FLT_MAX (updating minVal) or larger than -FLT_MAX (updating maxVal).

            - "minVal" = initialized to FLT_MAX, which is the largest possible value that can be represented by a float. Any number in the matrix will be smaller than this value, so minVal will be updated to the first value in the matrix. 
            - "maxVal" = initialized to -FLT_MAX, the smallest possible value. Any number in the matrix will be larger than this value, so maxVal will be updated to the first value in the matrix.
            
        Iterate through the 2D array:
            - The outer loop iterates over each row in the matrix
            - The inner loop iterates over each value in the row
            - For each value, code checks if its smaller than minVal or larger than maxVal (If either condition == true,  minVal or maxVal are updated accordingly.
        Output:
            After iteration through all values in matrix:
                - minVal contains the smallest value, and maxVal contains the largest value in the matrix.
        

        This technique is commonly used when you want to scan through a collection and find the min and max values without making any assumptions about the values in the collection.



Functions:

    readTiffImage():
        This function opens a TIFF file, retrieves image dimensions (width and height), as well as the total number of samples per pixel (samples_per_pixel), and reads the pixel data into a buffer (a 1D vector). Values are unsigned 8-bit integers from grayscale images (0-255).  Since my project focuses on implimenting my own image compression, encoding and decoding, this part of the code uses built-in functions from the libTIFF library (https://libtiff.gitlab.io/libtiff/) to open, and read image width, height, and pixel samples (samples per pixel).  
    
    Lempel-Ziv-Weltch (built-in algorithm for lossless tranforms):
        - Used in conjunction with other built-in read/write functions defined in the libTIFF library.  Essentially, after applying the wavelet transform and quantization, the resulting image still contains a lot of data. Even though the image is compressed by reducing the precision of the wavelet coefficients, there's still potential to make the file smaller.
    
    
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
The bitShiftAmount is used for applying bit-shifting quantization, which reduces the precision of image data by right-shifting the pixel values. This shift effectively narrows the range of pixel values, making it useful for compression by reducing the amount of data required to represent the image. 

    - bitShiftAmount = 0: No bit shift. Values are kept in full precision, and therefore no quantization is applied (image data remains.
    - bitShiftAmount = 1: Divides each pixel value by 2, reducing the range of pixel values.
    - bitShiftAmount = 2: Divides each pixel value by 4.
    - bitShiftAmount = 3: Divides each pixel value by 8, results in more aggressive compression.
    - bitShiftAmount = 4: Divides each pixel value by 16.
    - bitShiftAmount = 5: Divides each pixel value by 32.

As expected, higher bitShiftAmount values (e.g., 4 or 5) result in greater compression, but at the cost of reduced image quality due to more aggressive quantization. This can lead to a noticeable loss of detail, particularly in areas with fine textures or smooth gradients. On the other hand, lower bitShiftAmount values (e.g., 1 or 2) preserve more detail but achieve less compression. While this may not significantly reduce data storage, it still offers some level of quantization to decrease file size.  This can be demonstrated below:

    - bitShiftAmount = 1, or 2 (Light Compression: preserve image quality while still reducing the file size, useful for noise reduction)

![Image](images/v11_shift2_output_compressed.tif)

![Image](images/v11_shift2_output_recosntructed.tif)

    - bitShiftAmount = 3 (Moderate Compression: balances file size and image quality)

![Image](images/v11_shift3_output_compressed.tif)

![Image](images/v11_shift3_output_recosntructed.tif)

    - bitShiftAmount = 4, or 5 (Heavy Compression; reduce file size, more important than preserving image quality)
[TBD]

## Project revisions, and Future Development
Original thought: 
- Project data structure implemented:
    - Huffman tree or SPIHT tree (for encoding) (too much time considering the project scope)
    - bit shifting (for grayscale conversion) (wanted to use a method to convert both RGB images and grayscale images, and make use of bit shifting within this funcition, but went down rabbit hole of ARGB vs RGBA and also fell outside of main project scope)
        
[more needed, tbd]
