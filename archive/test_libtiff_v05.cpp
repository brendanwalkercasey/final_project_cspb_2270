#include <tiffio.h>
#include <iostream>
#include <vector>
#include <cstdint>  // For fixed-width integers
#include <cmath>
#define u32 std::uint32_t
#define u16 std::uint16_t
#define u8 std::uint8_t
using namespace std;

/*
v04 to v05 differences:
- removed function to convert RGB image to grayscale using the luminance formula
- added function for inverse transform of Mexican Hat Wavelet Transform:

       - Apply the inverse convolution: The inverse of a wavelet transform is typically done through a synthesis process, which is basically applying the inverse filter (a kernel) over the coefficients. In the case of the Mexican Hat wavelet, the inverse transformation involves reconvolving the coefficients with the wavelet kernel in reverse order. This step reconstructs the signal or image from the decomposed components.

       - Reconstruct the image: Once you've processed all coefficients through the inverse wavelet kernel, you can sum the results to get the reconstructed image or signal.

       - However, since Mexican Hat wavelets are non-orthogonal (they don't guarantee perfect reconstruction in the same way as orthogonal wavelets like Haar or Daubechies), achieving perfect reconstruction can be challenging and often depends on how the decomposition was performed.

- added "convertToUint8" = function is introduced to convert the double values (from the inverse transform) back to u8 values. This involves clamping the values to the valid 0-255 range, since pixel values should be integers between 0 and 255 for an 8-bit image.
    You can adjust the normalization or scaling if needed, depending on the characteristics of the transform.
- added reconstruction ima
*/

// Function to read the TIFF image into a buffer
// void readTiffImage(const char* filename, std::vector<std::uint8_t>& buffer, std::uint32_t& width, std::uint32_t& height, std::uint16_t& samples_per_pixel) {
void readTiffImage(const char* filename, vector<u8>& buffer, u32& width, u32& height, u16& samples_per_pixel) {      
    TIFF* tiff = TIFFOpen(filename, "r");
    if (tiff == nullptr) {
        cerr << "Could not open TIFF file!" << endl;
        return;
    }

    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);

    buffer.resize(width * height * samples_per_pixel);
    for (u32 row = 0; row < height; row++) {
        TIFFReadScanline(tiff, &buffer[row * width * samples_per_pixel], row);
    }

    TIFFClose(tiff);
}

// Function to compute the Mexican Hat (Ricker) wavelet at a point x
double mexicanHatWavelet(double x, double sigma) {
    return (1 - (x * x) / (sigma * sigma)) * exp(-x * x / (2 * sigma * sigma));
}


// Function to apply wavelet compression (simplified)
void waveletTransform(vector<u8>& buffer, u32 width, u32 height) {
    // Example: Apply a simple Haar wavelet transform (for illustration). 
    
    
    // Mexican Hat Wavelet Implimentation
    for (u32 row = 0; row < height; row++){
        for (u32 col = 0; col < width; col++){
            buffer[row * width + col] = mexicanHatWavelet(buffer[row * width + col], 10.0);
        }
    }

}

void inverseMexicanHatTransform(std::vector<double>& data, size_t width, size_t height, double sigma) {
    std::vector<double> result(data.size(), 0);

    // Mexican Hat kernel size, adjust based on your application
    int kernel_size = 7; // Example kernel size
    int offset = kernel_size / 2; // Kernel center

    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            double sum = 0;

            // Convolve with the Mexican Hat wavelet kernel
            for (int m = -offset; m <= offset; ++m) {
                for (int n = -offset; n <= offset; ++n) {
                    int x = j + n;
                    int y = i + m;

                    // Apply boundary checks
                    if (x >= 0 && x < width && y >= 0 && y < height) {
                        double kernel_value = mexicanHatWavelet(sqrt(m * m + n * n), sigma);
                        sum += data[y * width + x] * kernel_value;
                    }
                }
            }

            result[i * width + j] = sum;
        }
    }

    // Copy back the result to original data
    data = std::move(result);
}

// Function to convert double values to u8 (clamp to valid range)
void convertToUint8(const std::vector<double>& input, std::vector<u8>& output, size_t width, size_t height) {
    output.resize(width * height);

    for (size_t i = 0; i < width * height; ++i) {
        // Clamp the value to the range [0, 255]
        output[i] = static_cast<u8>(std::min(255.0, std::max(0.0, input[i])));
    }
}

// Function to write the TIFF image from the buffer
void writeTiffImage(const char* filename, const vector<u8>& buffer, u32 width, u32 height, u16 samples_per_pixel) {
    TIFF* tiff = TIFFOpen(filename, "w");
    if (tiff == nullptr) {
        cerr << "Could not open TIFF file for writing!" << endl;
        return;
    }

    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);  // Grayscale image will have 1 channel
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);  // 8 bits per pixel for grayscale
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);  // Grayscale, no RGB

    for (u32 row = 0; row < height; row++) {
        TIFFWriteScanline(tiff, const_cast<u8*>(&buffer[row * width]), row);
    }

    TIFFClose(tiff);
}

int main() {
    const char* inputFile = "board.tif";
    const char* outputFile = "output_decompressed.tif";

    u32 width, height;
    u16 samples_per_pixel;
    double sigma;
    vector<u8> buffer;

    // Step 1: Read the TIFF image
    readTiffImage(inputFile, buffer, width, height, samples_per_pixel);

    // Step 2: Apply wavelet compression (simplified)
    waveletTransform(buffer, width, height);

    // Step 3: Write the compressed image to a new TIFF file
    writeTiffImage(outputFile, buffer, width, height, samples_per_pixel);

    // Step 4: Inverse Wavelet Transform
    sigma = 1.0;

        // Perform the inverse transform
            // step 1: cast buffer to double vec
            vector<double> buff_double(buffer.begin(), buffer.end()); 
            // step 2: perfrom inverse transform
            inverseMexicanHatTransform(buff_double, width, height, sigma);

        //to see values, un-comment below:
        /*
        // Output the reconstructed image or signal
        for (size_t i = 0; i < height; ++i) {
            for (size_t j = 0; j < width; ++j) {
                std::cout << buff_double[i * width + j] << " ";
            }
            std::cout << std::endl;
        }
        */

    return 0;

    // Step 5: Convert the reconstructed data back to u8
    vector<u8> reconstructedImage;
    convertToUint8(buff_double, reconstructedImage, width, height);

    // Step 4: Write the reconstructed image to a new TIFF file
    writeTiffImage(outputFile, reconstructedImage, width, height, samples_per_pixel);

    cout << "Reconstructed image saved as " << outputFile << endl;

    cout << "Wavelet compression completed!" << endl;

    return 0;
}