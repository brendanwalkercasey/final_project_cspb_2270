#include <tiffio.h>
#include <iostream>
#include <vector>
#include <cstdint>  // For fixed-width integers
#include <cmath>
//#define M_PI 3.14159265358979323846
#define u32 std::uint32_t
#define u16 std::uint16_t
#define u8 std::uint8_t
using namespace std;

/*
Key changes from v06 and v06-1 (test version to "normalize values"):
- essentially a modified v03
- commented out rgbToGrayscale() method (already using grayscale image)
- Haar instead of Mexican Hat Wavelet Transform
- explicitly declared samples per pixel in main function
    - You were initially passing samples_per_pixel to the writeTiffImage function because the input TIFF image could be RGB or grayscale, and this value would change depending on the image type. When the image is converted to grayscale, the samples_per_pixel would be set to 1.
    - However, after converting to grayscale, you should still account for samples_per_pixel in the TIFF header, even if the image becomes a single-channel grayscale image. The samples_per_pixel value in the TIFF header indicates how many components (or channels) each pixel contains, and it needs to reflect the actual image data format.
- Try COMPRESSION_LZW: The LZW compression is a lossless compression algorithm that works well for images with large areas of uniform color.
- tried different image to see different results
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


// Function to convert RGB image to grayscale using the luminance formula
void rgbToGrayscale(std::vector<u8>& buffer, u32 width, u32 height, u32 samples_per_pixel) {
    // Make a new buffer to store the grayscale image
    std::vector<u8> grayscaleBuffer(width * height);

    for (u32 i = 0; i < height; i++) {
        for (u32 j = 0; j < width; j++) {
            // Extract the RGB pixel (Assuming 8-bit per channel RGB image)
            u32 pixel = buffer[i * width * samples_per_pixel + j * samples_per_pixel];

            // Extract RGB components
            u8 red = (pixel >> 16) & 0xFF;
            u8 green = (pixel >> 8) & 0xFF;
            u8 blue = pixel & 0xFF;

            // Compute the grayscale value using the luminance formula
            u8 grayscale = static_cast<u8>(red * 0.299 + green * 0.587 + blue * 0.114);

            // Store the grayscale value in the new buffer (one value per pixel)
            grayscaleBuffer[i * width + j] = grayscale;
        }
    }

    // Copy the grayscale buffer back to the original buffer (if you want to overwrite the original image)
    buffer = std::move(grayscaleBuffer);
}



// Mexican Hat Wavelet Transform
double mexicanHatWavelet(double num, double sigma) {
    return (2.0 / (sqrt(3) * pow(M_PI, 0.25))) * (1 - num * num / (sigma * sigma)) * exp(-num * num / (2.0 * sigma * sigma));
}

// Function to apply inverse Mexican Hat wavelet transform (decompression)
void inverseMexicanHatTransform(vector<double>& data, size_t width, size_t height, double sigma) {
    vector<double> result(data.size(), 0);
    int kernel_size = 7; // Example kernel size
    int offset = kernel_size / 2; // Kernel center

    // Apply the inverse transform with the Mexican Hat wavelet
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            double sum = 0;
            for (int m = -offset; m <= offset; ++m) {
                for (int n = -offset; n <= offset; ++n) {
                    int x = j + n;
                    int y = i + m;
                    if (x >= 0 && x < width && y >= 0 && y < height) {
                        double kernel_value = mexicanHatWavelet(sqrt(m * m + n * n), sigma);
                        sum += data[y * width + x] * kernel_value;
                    }
                }
            }
            result[i * width + j] = sum;
        }
    }
    data = std::move(result); // Update the data with the decompressed result
}


// Function to apply wavelet compression (simplified)
// We modify the buffer to use double precision for wavelet operations,
// and normalize pixel values to fall in 8-bit range
void waveletTransform(vector<u8>& buffer, u32 width, u32 height, double sigma) {

    // Create a double precision buffer for wavelet transform
    vector<double> tempBuffer(width * height);

    // Convert the u8 buffer to double for the wavelet transform
    for (u32 i = 0; i < width * height; i++) {
        tempBuffer[i] = static_cast<double>(buffer[i]);
    }

    // Apply the Mexican Hat Wavelet (using tempBuffer for double precision)
    for (u32 row = 0; row < height; row++) {
        for (u32 col = 0; col < width; col++) {
            double pixelValue = tempBuffer[row * width + col];
            tempBuffer[row * width + col] = mexicanHatWavelet(pixelValue, sigma); // Pass sigma here
        }
    }

    // Normalize the result to [0, 255] range for output (8-bit)
    double minVal = *std::min_element(tempBuffer.begin(), tempBuffer.end());
    double maxVal = *std::max_element(tempBuffer.begin(), tempBuffer.end());
    
    for (u32 i = 0; i < width * height; i++) {
        // Normalize to 0-255 range
        double normalized = (tempBuffer[i] - minVal) / (maxVal - minVal) * 255.0;
        buffer[i] = static_cast<u8>(std::clamp(normalized, 0.0, 255.0));
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

    // Set the compression method to lossless (e.g., LZW or Deflate)
    TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_LZW);  // You can also try COMPRESSION_DEFLATE or COMPRESSION_PACKBITS

    for (u32 row = 0; row < height; row++) {
        TIFFWriteScanline(tiff, const_cast<u8*>(&buffer[row * width]), row);
    }

    TIFFClose(tiff);
}


int main() {
    const char* inputFile = "cameraman.tif";
    const char* outputFile = "output_compressed.tif";

    u32 width, height;
    u16 samples_per_pixel = 1;
    double sigma = 99; // Example sigma, change as needed
    vector<u8> buffer;

    // Step 1: Read the TIFF image
    readTiffImage(inputFile, buffer, width, height, samples_per_pixel);

    // Step 2: Convert RGB to Grayscale (if applicable)
    rgbToGrayscale(buffer, width, height, samples_per_pixel);

    // Step 3: Apply wavelet compression (simplified)
    waveletTransform(buffer, width, height, sigma);

    // Step 4: Write the compressed grayscale image to a new TIFF file
    writeTiffImage(outputFile, buffer, width, height, samples_per_pixel);
    
    // Step 5: Apply the inverse transform to reconstruct the image (optional)
    vector<double> buff_double(width * height); // Use a double buffer for inverse transform
    inverseMexicanHatTransform(buff_double, width, height, sigma); // Use an appropriate sigma value

    cout << "Wavelet compression completed!" << endl;

    return 0;
}
