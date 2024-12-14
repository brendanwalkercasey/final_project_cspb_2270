# Wavelet Compression of Imagery
Final Project from CU Boulder's Data Structures 2270, under Applied Computer Science Bachelors Program.


## Overview:
Wavelet compression is a popular method for image compression, which leverages wavelet transforms to represent an image in a more efficient way, typically achieving high compression ratios with minimal loss of quality. The general idea is to decompose the image into a set of wavelet coefficients that can be quantized and encoded more effectively than the original pixel values.  The wavelet transform is a mathematical tool that decomposes an image into different frequency components at multiple scales. These components capture both the low-frequency information (which contains the smooth, overall features of the image) and high-frequency details (which capture sharp edges and fine textures).

In simple terms, it’s akin to applying a filter to the image at different scales, allowing you to see both coarse and fine details in the image. By representing the image in terms of wavelet coefficients, we can selectively retain the most significant information while discarding less important data.

This project uses a Daubechies Wavlet (Daubechies-4): for decomposition of an image into frequency components, bit-shifting quantization: to quantize large values, while smaller coefficients are thresholded out with a sigma value, an inverse wavelet transform: to reconstruct the image from its compressed components, and built-in methods (from libTIFF library) to encode the image after quantization (with Limpel-Ziv-Welch, or LZW Compression).  

## Data Structures/Techniques Implemented:

    - Vectors (main base structure for storing image data, in the form of 1-D and 2-D arrays, typically with 8-bit, 16-bit, unsigned integers, or float types)

    - Buffers (temporary structures used during 2D Daubechies tranform to store intermediate results)

    - Pointers (namely unique pointers to TIFF objects, to interact with the libTIFF library in a memory efficent and memory safe manner)

    - Standard Library Algorithms (standard built-in methods to better organize code, and save time e.g. min_element and max_element for finding min/max values)

    - Bit-Shifting (simple and fast operations with minimal computations to reduce data on individual pixel values, effectively compressing image file with fewer bits)



## Variables:

    Wavelet Decomposition (low/high pass filters): The Daubechies-4 filter is used for the 2D wavelet transform. The decomposition into sub-bands happens in two steps: applying the wavelet transform to the rows and then to the columns:

        - Low-pass filter (h0, h1, h2, h3)
        - High-pass filter (g0, g1, g2, g3)
        - Inverse filter parmeters (low) (ih0, ih1, ih2, ih3) 
        - Inverse filter parmeters (high) (ig0, ig1, ig2, ig3)

    These are the inverse Daubechies-4 wavelet transform coefficients (used for reconstructing the image after compression).
        
        - For each row, the daubechies1D function is called, which applies the Daubechies-4 filter (low-pass and high-pass) to the data. After processing, the data is split into two parts: 1) low frequencies and 2) high frequencies.
        The same process happens for every column in the daubechies2D function, further splitting the image into frequency bands.

    Input Parameters:
        - inputFile: string representing the path to the input TIFF file to be read ("cameraman.tif").
        - outputFile: string representing the path where the compressed image will be saved.
        - outputReconstructedFile: string representing the path where the reconstructed image after applying the inverse wavelet transform will be saved.
        - sigma: double representing the threshold for soft thresholding in the quantization step. This controls how aggressively small values in the transformed image are set to zero during quantization.
        - bitShiftAmount: integer representing the amount of bit shifting to apply during quantization (e.g., shifting by 2 bits). This reduces the precision of the wavelet coefficients, simulating compression.
    
    Image and Buffer Information:
        - width: uint32_t variable representing the width of the image (in pixels). This is determined by reading the image from the TIFF file.
        - height: uint32_t variable representing the height of the image (in pixels). This is also determined by reading the image from the TIFF file.
        - samples_per_pixel: uint16_t variable representing the number of color channels per pixel. For grayscale images, this is typically 1 (single channel).
        - buffer: std::vector<uint8_t> (vector of unsigned 8-bit integers) used to store the pixel data read from the TIFF file or written to the TIFF file. This holds image data (grayscale) in a 1D array format after reading the TIFF file. During processing, it's resized and updated with  transformed, quantized, or reconstructed pixel values.
    
    Temporary Buffers and Containers:
        - tempBuffer: vector<double> used as temporary buffer to hold the image data as double values during the wavelet transform and reconstruction processes. This allows for more precise computations before converting the data back to uint8_t for output.

    Functions and Flow Control
        - TIFF* tiff: pointer to a TIFF structure used in conjunction with the TIFFOpen and TIFFClose functions for reading and writing TIFF files. The std::unique_ptr ensures that the TIFF file is properly closed when it goes out of scope.
        - temp: Used for daubechies1D() and inverseDaubechies1D() functions to store temporary transformed data during the forward and inverse Daubechies-4 transform steps.
        - normalized: double used for normalizing the pixel values (for example, during the inverse wavelet transform and when preparing the buffer for output), ensuring they fit in the 0-255 range for proper display in an 8-bit grayscale image.


    Sub-Band Parameters: 

        minVal and maxVal: The main reason for using constants is for defining initial extremes for the minVal and maxVal variables. When starting with largest possible, and smallest possible float values, we ensure that every value in the collection will either be smaller than FLT_MAX (updating minVal) or larger than -FLT_MAX (updating maxVal).

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



## Functions:

    readTiffImage():
        - Opens TIFF image file for reading and loads the image data into a buffer. It extracts key metadata such as the image width, height, and the number of samples per pixel (such as grayscale or color channels). After resizing the buffer to fit the image dimensions, it reads the image scanlines (rows of pixels) one by one and stores them in the buffer. If the file cannot be opened, an error message is displayed. Values are unsigned 8-bit integers from grayscale images (0-255).  Since my project focuses on implimenting my own image compression, encoding and decoding, this part of the code uses built-in functions from the libTIFF library (https://libtiff.gitlab.io/libtiff/) to open, and read image width, height, and pixel samples (samples per pixel).  

    daubechies1D():
        - Performs 1D discrete wavelet transform (DWT) using Daubechies-4 coefficients on a vector of data. It applies both low-pass and high-pass filters to consecutive data values, separating the signal into approximation (low frequency) and detail (high frequency) components. These components are then stored in a temporary vector (temp) and used to update the original data, which is transformed in-place. 
    
    daubechies2D():
        - Applies 1D Daubechies-4 wavelet transform (calculated previously) to a 2D image, transforming each row of the image, then applying the same transformation to each column. The process effectively decomposes the image into various frequency components, capturing both horizontal and vertical details. This function modifies the image in-place to store the transformed data.

    applyBitShiftingQuantization():
        - Quantizes data using a simple bit-shifting technique. It iterates through each value in the data and applies soft thresholding by setting values below a given threshold (sigma) to zero. For values above the threshold, the function applies a bit-shift (by 2 bits in this case) to reduce the precision, effectively compressing the data (In essence, sigma controls which values are kept (by zeroing out small ones), while bitShiftAmount controls how precisely the kept values are represented (by reducing their bit precision). Both together help in compressing the image data). The result is a quantized version of the original data with reduced size.  For the sake of simplicity, I decided to shift by 2 bits.

    waveletTransform():
        - Applies the Daubechies wavelet transform to decompose the image into multiple subbands (such as approximation, horizontal details, vertical details, and diagonal details). The basic idea behind applying Daubechies wavelets for image compression is to use the wavelet transform to decompose the image into different frequency subbands. These subbands are then quantized and encoded, with the most significant subbands preserved to maintain image quality. The less significant parts (those representing high-frequency noise or detail) can be discarded to achieve compression.  Daubechies wavelets are a family of orthogonal wavelets (meaning that the wavelet functions are mutually orthogonal), developed by Ingrid Daubechies in the late 1980s. They are widely used in signal processing and image compression because of their compact support (non-zero values only within a limited region) and smoothness. These properties make them ideal for efficiently representing signals (especially images) while maintaining a good balance between compression and reconstruction quality. Daubechies wavelets are defined by the number of vanishing moments they have, which influences their ability to approximate polynomial functions. For example, DbN wavelets (Daubechies wavelet with N vanishing moments) are characterized by the number of moments the wavelet function's integral is zero. Db4, Db6, Db8, etc., are common examples, where the number indicates the number of vanishing moments. My code uses 4 total.
    
    writeTiffImage():
        - Writes an image buffer to a TIFF file, setting various properties of the image such as width, height, and number of samples per pixel. The image data is written row-by-row as scanlines, and compression (LZW) is applied to reduce file size. This function also ensures that the image is saved in a grayscale format with 8 bits per sample.

    inverseDaubechies1D(): 
        - The inverseDaubechies1D function reverses the effect of the Daubechies-4 1D wavelet transform. It reconstructs the original signal by applying inverse low-pass and high-pass filters to the transformed data. The function operates on the transformed data in-place, updating it with the inverse transformation result, effectively reversing the approximation and detail separation done by the forward wavelet transform.

    inverseDaubechies2D(): 
        - Applies inverse of the 1D Daubechies-4 wavelet transform to both the rows and columns of a 2D image matrix. It reconstructs the original image by applying the inverse wavelet transform in a manner similar to daubechies2D, first to each row and then to each column, effectively merging the approximation and detail components. The result is an image that approximates the original before any transformation was applied. [Need more...]

    inverseWaveletTransform(): 
        - Applies inverse wavelet transform to an image buffer that was previously compressed using the Daubechies-4 wavelet transform. The image is first converted to a temporary double-precision buffer, then the inverse 2D Daubechies-4 transform is applied. The image is then normalized back to the 8-bit range and written to the buffer, reconstructing the image as closely as possible to its original form before compression.

    writeSubbandToTiff(): 
        - Writes specific sub-band of a wavelet-transformed image (such as LL, LH, HL, or HH) to a TIFF file. It first finds the minimum and maximum values in the sub-band, normalizes the data to an 8-bit grayscale range (0–255), and then writes sub-band to a TIFF file. This function allows for visualizing individual wavelet sub-bands, which contain different frequency components of the image.

    daubechies2DWithSubbands(): 
        - Decomposes a 2D image into four distinct sub-bands (LL, LH, HL, HH) using the Daubechies-4 wavelet transform. It applies the 1D transform to each row and column of the image, and then separates the result into approximation (LL) and detail (LH, HL, HH) components. These sub-bands represent different frequency content of the image: low frequencies (LL), horizontal details (LH), vertical details (HL), and diagonal details (HH). [This decomposition is useful for image compression and analysis...need more]

    Lempel-Ziv-Weltch() (built-in algorithm for lossless tranforms):
        - Used in conjunction with other built-in read/write functions defined in the libTIFF library.  Essentially, after applying the wavelet transform and quantization, the resulting image still contains a lot of data. Even though the image is compressed by reducing the precision of the wavelet coefficients, there's still potential to make the file smaller.

## Steps:

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
        
        Additional Analysis: Using sigma as a threshold to remove smaller wavelet coefficients before quantization is useful for reducing noise or compressing the image more aggressively.  I decided to add a thresholding operation based on sigma before applying the quantization. Values smaller than sigma will be set to zero, which could be useful if sigma represents a specific noise level we can identify for suppression.  

    4) Reconstruction (Inverse Wavelet Transform):

        Action: The inverse Daubechies-4 wavelet transform is applied to reconstruct the image from its wavelet sub-bands.

        Details: The image is reconstructed by applying the inverse wavelet transform (first on columns, then on rows), reversing the effect of the initial wavelet transform.

    5) Post-process and Encode:

        5a) Normalize the Data: The transformed image data is normalized back to the 0-255 grayscale range for proper display. This normalization ensures that the pixel values fit within the standard 8-bit range.

        5b) Save the Processed Image: The compressed (and potentially quantized) image is saved to a new TIFF file using LZW compression for file size reduction.
        
        5c) Reconstruct and Save the Image: The image is reconstructed using the inverse wavelet transform (as described earlier) and then saved to another TIFF file as the "reconstructed" version.

    6) Additional Steps: 

        - Sub-Band Decomposition: After applying the wavelet transform, four sub-bands are created to illustrate details of directional components of the image:

            - LL (approximation) (visible in top left of compressed output tiff)

            - LH (horizontal details) (visible in bottom left of compressed output tiff)

            - HL (vertical details) (visible in top right of compressed output tiff)

            - HH (diagonal details) (visible in bottom right of compressed output tiff)

        While these details are visible in the compressed file (see list above), the code for generating individual output tiff files (e.g. LL tiff, or HH tiff) did not correctly split and resize the original compressed output image.  Therefore, generating individual tiff files for each sub-band was left as optional under the main() function.  


    
## How to run code:

    required installations:
        - libTIFF library (https://libtiff.gitlab.io/libtiff/)
        - GNU Compiler Collection (https://gcc.gnu.org/)

    - My code can be compiled and run as follows:
        
        Compile Command (with flags):

            g++ -std=c++17 -o wavelet_compress wavelet_compression.cpp -ltiff

        Run Executable:

            ./wavelet_compress

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


## Analysis:

Bit Shift Values: The bitShiftAmount is used for applying bit-shifting quantization, which reduces the precision of image data by right-shifting the pixel values. This shift effectively narrows the range of pixel values, making it useful for compression by reducing the amount of data required to represent the image. 

    - bitShiftAmount = 0: No bit shift. (values are kept in full precision, and therefore no quantization is applied, and image data remains.)
    - bitShiftAmount = 1: Divides each pixel value by 2.
    - bitShiftAmount = 2: Divides each pixel value by 4. (reduces the range of pixel values.)
    - bitShiftAmount = 3: Divides each pixel value by 8.
    - bitShiftAmount = 4: Divides each pixel value by 16. (results in more aggressive compression.)
    - bitShiftAmount = 5: Divides each pixel value by 32.
    - bitShiftAmount = 6: Divides each pixel value by 64. (greatly reduces image quality and size, but preserves most basic features.)

As expected, higher bitShiftAmount values (e.g., 4 or 5) result in greater compression, but at the cost of reduced image quality due to more aggressive quantization. This can lead to a noticeable loss of detail, particularly in areas with fine textures or smooth gradients. On the other hand, lower bitShiftAmount values (e.g., 1 or 2) preserve more detail but achieve less compression. While this may not significantly reduce data storage, it still offers some level of quantization to decrease file size.  The most interesting was the visual result of each of the decomposed sections (i.e. sub-bands: HL, LH, HH, and LL) of the wavelet_compress_shift6_sigm42_reconstructed.tif, when increasing the bit shifting value from 4 to 6.  The decomposed sections of horizontal, vertical and diagonal components show stronger edges than previously demonstrated (i.e. stronger edges when shifting bits by 6 as opposed to 4 or 2).  Additionally, there is a significant amount detail still preserved when compiling and creating the wavelet_compress_shift6_sigm42_reconstructed.tif.  This can be demonstrated below:

    - bitShiftAmount = 1, or 2 (Light Compression: preserve image quality while still reducing the file size, useful for noise reduction)

![Image](outputs/wavelet_compress_shift2_sigm42_compressed.tif) (29KB)

![Image](outputs/wavelet_compress_shift2_sigm42_reconstructed.tif) (38KB)

    - bitShiftAmount = 3 (Moderate Compression: balances file size and image quality)

![Image](outputs/wavelet_compress_shift4_sigm42_compressed.tif) (15KB)

![Image](outputs/wavelet_compress_shift4_sigm42_reconstructed.tif) (31KB)

    - bitShiftAmount = 4, 5, or 6 (Heavy Compression; significantly reduced file size, where size reduction is more important than preserving image quality, and edges more pronounced in individual decomposed horizontal, vertical and diagonal components)

![Image](outputs/wavewavelet_compress_shift6_sigm42_compressed.tif) (7KB)

![Image](outputs/wavelet_compress_shift6_sigm42_reconstructed.tif) (18KB)


Daubechies-4 Algorithm: As expected, the Daubechies Compression breaks down the image into "approximation" and "detail" coefficients (low and high-frequency components), which helps represent the image (e.g. cameraman.tif) more compactly by separating out important structures (like hard edges and textures of the man's sillouete, his tripod, camera, and background building features) from less important ones (like smooth regions of the sky, grass, and repeating patterns of vegetation).  While I experiemented with vastly different sigma values, this didnt seem to change the image reconstruction from one iteration to the next, and unexpectedly, sizes were virtually the same between significantly different sigma values (e.g. wavelet_compress_shift2_sigma3_reconstructed.tif and wavelet_compress_shift2_sigma99_reconstructed.tif were both 38KB).  This is likely due to improper implementation of my sigma value, or perhaps improper combination with other compression techniques used (combined with bitshifing and/or my implimentation of the Lempel-Ziv-Weltch() function).  However this implimentation was successful in de-noising the image from specific artifacts (e.g. removing the random white band artifact in the bottom half of the cameraman.tif, and reducing complexity of pixel values in the sky in the grass). 

## Future Development
    Ongoing Issues:
    1) Intermittent Data Loss in Output Files: The largest issue this code faces is an unknown, artifact that causes an occasional memory leak.  Admittedly this was difficult to diagnose, and given more time, I plan to find and fix this issue.  Essentially, about 1 and every 5 iterations of compiling the the code and running the "wavelet_compress" binary, it produces two (blank) grey or black imagery with one or two white bands (narrow and horizontal), instead of a succuessfully uncompressed image and a 4-band decomposition image.  I spent sevearal days debugging my code to add logic that would (hopefully) combat this issue.  These efforts included:
            - re-factoring my pointer structure (adding custom deleter fucntion for TIFF* to ensure proper cleanup, and using std::unique_ptr with custom deleter to automatically close the TIFF file)
            - Ensuring image width and height are even for sub-band splitting (both in readTiffImage() and waveletTransform()), effectively creating a "padded" width/height by rounding up to next multiple of 4
            - adjusting bit shift value in applyBitShiftingQuantization(), since higher bit shift values may lead to significant data loss.
            - additional images were used (of varying sizes), but the issue persisted, leading me to beleive the bug still lies with how the code impliments the transform, resizing, or compression.
            - explicitly resizing the buffer with padded width and height in my waveletTransform(), and padding image with zeros.
            - Changing the "#define" to "typedef" when aliasing my uint data types (e.g. uint8_t as u8), initially beleiving it was an issue at the preprocessor level, affecting the entire translation unit (file).  Then removing (commenting out) typedef entirely, and falling back on using regular data types.  The initial thought behind aliasing my data types were simply to improve readability of my code.

        Ultimatley, I beleive the issue is some type of memory leak, or incorrect/overcompression technique, which results in blank grey images periodically in my code's output (depending on the image input).  Perhaps the combination of bit shifting and wavelet compression requires special sigma values, or customized Lempel-Ziv-Weltch() function to be calculated depending on image input and its spatial patterns.  Currently, in applyBitShiftingQuantization(), I apply bit-shifting to compress the image. The current shift is 2 bits, which may lead to significant data loss. While I've adjusted this value based on the level of compression I may want to also consider other quantization techniques and/or use more sophisticated approaches.  All efforts above are documented in "wavelet_compression_debug_grey.cpp".

    2) Improperly Sized Sub-Band Output Files:  While the above issue took ~90% of my time allocated towards debugging, I also experienced issues with the sub-band outputs incorrectly displaying individual decomposed sub-bands for Vertical, Horizontal, Diagonal, and Approximation components (i.e. sub-bands: HL, LH, HH, and LL).  This is an issue that I believe is an easier fix, however my priority was diagnosing an adding logic to debug the main file outputs (reconstructed image after transform and compression, and decomposed output image with all four sub-bands).

## Original Plans: 
    - OpenCV Library vs. LibTIFF Libraries: My original plan was to use OpenCV for its ease of use (and use LibTIFF as a backup library).  However, the route I took to install openCV (both using package managers, and precompiled binaries for quick installs) required several dependencies:

                - Package managers:
                    - Homebrew package manager (decided not to use, until I had more information about this issue - see https://saagarjha.com/blog/2019/04/26/thoughts-on-macos-package-managers/)
                    - Macports package manager (decided not to use - see above, requires a lot of manual configurations, which runned the risk of time spent learning macports vs data structures needed for my project)

                - Examples of General Dependencies from downloading from source:
                    - pkg-config 
                        - when i ran pkg-config --modversion opencv4, my terminal returned zsh: command not found: pkg-config
                    - pkg-config required glib which depending on the version downloaded, could not run ./confgure, and instead required meson tool

                - Issues downloading from pre-compiled binaries:
                    - opencv pre-compiled binaries wouldn't link to my test_opencv.cpp code (test code for reading an image)
        
        - Libraries Summary: Logistically, for the scope of my project, it made more sense to use libTIFF library over OpenCV's library.  From my research, OpenCV is designed for software developers who dont want to hassle with the complexities of low-level file handling, or who are working with several image formats, and don’t need to focus exclusively on TIFF files. Using libTIFF turned out to be beneficial because it is designed for TIFFs specifically, and the API for libtiff is lower-level/ more complex than OpenCV’s. Code built on libTIFF source code needs to manage memory allocation, handle different compression types, and deal with the intricacies of the TIFF format manually.  This library provided a more in depth experience for a final project in Data Structures because my project requires fine-grained control over TIFF files, and specific compression methods.

    - Additional Data Structure(s): 
    
        - Bit Shifting for Grayscale Conversion: I wanted to use a method to convert both RGB images and grayscale images, and make use of bit shifting within this function.  I went down rabbit hole of ARGB vs RGBA methods, and how bit shifting and bit masking are applied in different contexts, especially when dealing with different pixel formats and how the individual color channels (Red, Green, Blue, Alpha) are stored in memory.  While it was benefical to learn about RGBA (0xAARRGGBB) (where the red channel is in bits 16–23, so you shift 16 bits to the right to extract it) as well as ARGB (0xRRGGBBAA) (where the red channel is in bits 24–31, so you shift 24 bits to the right to extract it), I determined this fell outside of main project scope, where this time could be better used to learn about wavelet transforms, and apply bitshifting towards compression techniques (the main focus of my project).  

        - SPIHT (Set Partitioning in Hierarchical Trees) is an algorithm specifically designed for wavelet-based image compression, and it's a popular choice because it achieves high compression ratios while preserving image quality. It operates by progressively encoding the most significant wavelet coefficients first, and is particularly well-suited for progressive compression, where an image can be progressively reconstructed as more bits are received.  This was an interesting direction, but was ultimatley abandoned due to time constraints.  However implementing SPIHT fully in C++ requires more attention to data structures like trees and priority queues, and efficient handling of the bitstream for encoding, which would be a great follow-on to this project.
    
    - Additional Data Sets:
        
        - In my proposal, I planned to compare a wavelet transform effect on medical imagery with landsat data (individual RGB bands for simplicity), however considering the scope of this project, I re-evaluated chose to focus on implimenting the compression algorithms (specifically learning the Daubechies-4 method, as proposed) to incorperate the necessary (and desired) data structures (e.g. bit-shifting, low-pass filters, re-constructing the image, etc.).  Therefore I decided to use a single TIFF image and break down its components into the frequencies of its pixel values, and illustrate directional sub-bands.  Additionally, if my objective were to analyze the content of the Landsat image (e.g., vegetation analysis, land use classification, etc.), it was important to consider the raw compression techniques (wavelet transform + SPIHT) might not be the best option unless I was focusing on preserving high-level features. Lossy compression with no additional refinemnt (e.g. atmospheric correction, etc.) can introduce artifacts, which could impact feature extraction accuracy.


## Notes:

    - Additional Compilation Errors:
        - #define M_PI 3.14159265358979323846 commented out because it caught its orginal definition when compiling (rendering this uneccessary).

## Resources

1. TIFF Manipulation (I/O, and built in image processing functions):

    - LibTIFF library: https://wavelet2d.sourceforge.net/
    - LibTIFF Manual: http://www.libtiff.org/man.html
    - https://pywavelets.readthedocs.io/en/latest/ref/wavelets.html
    - https://graphics.stanford.edu/wikis/cs148-07/Assignment7
    - Texas A&M lecture: https://people.qatar.tamu.edu/tingwen.huang/math414/5.1.pdf
    - Open Source tiff files: https://people.math.sc.edu/Burkardt/data/tif/tif.html

2. Daubechies Wavelets and Wavelet Transforms

    - A Primer on Wavelets and Their Scientific Applications by James S. Walker (Portion of Book): https://books.google.com/books?id=F-9VGmsdbdYC&pg=PA29&source=gbs_toc_r&cad=2#v=onepage&q&f=false
    - http://bearcave.com/misl/misl_tech/wavelets/index.html
    - http://bearcave.com/misl/misl_tech/wavelets/daubechies/index.html
    - http://bearcave.com/software/java/wavelets/daubechies/index.html
    - http://bearcave.com/misl/misl_tech/wavelets/daubechies/daub.h


3. Wavelet Compression Theory

    - https://kids.frontiersin.org/articles/10.3389/frym.2023.1200611
    - https://www.nytimes.com/2021/09/14/magazine/ingrid-daubechies.html
    - https://www.researchgate.net/profile/Don-Hong/publication/227246220_Wavelet_Image_Compressor-MinImage/links/64fb85fb90dfd95af61ff2a4/Wavelet-Image-Compressor-MinImage.pdf
    - Grgic, S., Kers, K., & Grgic, M. (1999, July). Image compression using wavelets. In ISIE'99. Proceedings of the IEEE International Symposium on Industrial Electronics (Cat. No. 99TH8465) (Vol. 1, pp. 99-104). IEEE.
    - Adriaens, J., & Palsetia, D. SIMD Implementation of the Discrete Wavelet Transform.
    - Mohammed, F. G., & Al-Dabbas, H. M. (2018). The Effect of Wavelet Coefficient Reduction on Image Compression Using DWT and Daubechies Wavelet Transform. Science International, 30(5), 757-762.
    - A Review of Image Compression Techniques by Rajandeep Kaur and Pooja Choudhary (Research Paper): https://d1wqtxts1xzle7.cloudfront.net/46755270/rajan_pooja-libre.pdf?1466766740=&response-content-disposition=inline%3B+filename%3DA_Review_of_Image_Compression_Techniques.pdf&Expires=1733956989&Signature=Mna-hra46L4bAr9ET8SXYAQ-EaFNhmGKHx1OCSSv4UmIQ-qVedrQM04dNo304HV5qRSfsJe7Yn~B5617h0ogcliYphOiKJL8vO0R5JOpLZQxtHP08Fsni6-uHb9uhREhGzwTiE3wOL7-7LwDBs4Loso8BM8~xtTX8p~f6cKjHYy5~~gmmbeQS8TpmzmZtSklEwvab73D-By1ZSokL~-AtOs5anOyaRqfbqjsOeped3JKZgkipkJSuY6HV-cPmF72iAG5447B0WJnLPMIxM2plgkm76fOo1o6j76qPZuNJz1bzUTVaCGp-P5Za7yvp7rS3~n7AHL4hrHkDZMHsgJ5cA__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA 
    - https://dc.etsu.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=1108&context=etd
    - http://www.stat.columbia.edu/~jakulin/Wavelets/index.html
    - https://ntrs.nasa.gov/api/citations/19660001929/downloads/19660001929.pdf
    - https://iopscience.iop.org/article/10.1088/1755-1315/280/1/012031/pdf

4. Bit Shifting Quantization:
    
    - https://www.learncpp.com/cpp-tutorial/bit-manipulation-with-bitwise-operators-and-bit-masks/
    - Bitwise Operators in C++ (Tutorial): https://www.learncpp.com/cpp-tutorial/bitwise-operators/
    - C++ Bitwise Operators Cheat Sheet: https://en.cppreference.com/w/cpp/language/operator_arithmetic

5. General Image Processing Resources:

    - OpenCV Documentation: https://docs.opencv.org/master/ (ended up not using, but referenced for content and image processing concepts)
    - https://stackoverflow.com/questions/60377598/change-rgb-to-grayscale-in-c (ended up ot using in final iteration, and just used grayscale images, but referenced for content and concepts for bitshifting)
    - https://it.nc.gov/documents/files/understanding-compression-geospatial-raster-imagery/open (source for idea to use Lempel-Ziv-Welsh (LZW) method)

6. Visualization of Wavelet Transforms

    - Matplotlib for Wavelets (https://www.mathworks.com/help/wavelet/gs/introduction-to-the-wavelet-families.html)