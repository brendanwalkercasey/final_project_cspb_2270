#include <tiffio.h>
#include <iostream>
#include <vector>
#include <cstdint>  // For fixed-width integers
#include <cmath>
#include <algorithm> // For std::min, std::max
#define u32 std::uint32_t
#define u16 std::uint16_t
#define u8 std::uint8_t
using namespace std;


/*
Key changes from v05:
- Perform the inverse wavelet transform on the image after compression.
- Convert the resulting double values back to u8 for proper image storage.
- Save the decompressed image to a new TIFF file.
*/

// Function to read a TIFF image
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
    // Apply a simple Mexican Hat wavelet transform (this is for compression simulation)
    for (u32 row = 0; row < height; row++) {
        for (u32 col = 0; col < width; col++) {
            buffer[row * width + col] = static_cast<u8>(mexicanHatWavelet(buffer[row * width + col], 10.0));  // Simplified transform
        }
    }
}

// Function to apply inverse Mexican Hat wavelet transform (decompression)
void inverseMexicanHatTransform(vector<double>& data, size_t width, size_t height, double sigma) {
    vector<double> result(data.size(), 0);
    int kernel_size = 7; // Example kernel size
    int offset = kernel_size / 2; // Kernel center

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

// Function to convert double values to u8 (clamp to valid range)
void convertToUint8(const vector<double>& input, vector<u8>& output, size_t width, size_t height) {
    output.resize(width * height);
    for (size_t i = 0; i < width * height; ++i) {
        output[i] = static_cast<u8>(std::min(255.0, std::max(0.0, input[i])));
    }
}

// Function to write a TIFF image from the buffer
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
    const char* compressedFile = "output_compressed.tif";

    u32 width, height;
    u16 samples_per_pixel;
    double sigma = 1.0;
    vector<u8> buffer;

    // Step 1: Read the TIFF image
    readTiffImage(inputFile, buffer, width, height, samples_per_pixel);

    // Step 2: Apply wavelet compression (simplified)
    waveletTransform(buffer, width, height);

    // Step 3: Write the compressed image to a new TIFF file
    writeTiffImage(compressedFile, buffer, width, height, samples_per_pixel);

    // Step 4: Inverse Wavelet Transform
    // Convert the buffer to double type for the inverse transform
    vector<double> buff_double(buffer.begin(), buffer.end());

    // Perform the inverse transform
    inverseMexicanHatTransform(buff_double, width, height, sigma);

    // Step 5: Convert the reconstructed data back to u8
    vector<u8> reconstructedImage;
    convertToUint8(buff_double, reconstructedImage, width, height);

    // Step 6: Write the decompressed image to a new TIFF file
    writeTiffImage(outputFile, reconstructedImage, width, height, samples_per_pixel);

    cout << "Reconstructed image saved as " << outputFile << endl;
    cout << "Wavelet compression and decompression completed!" << endl;

    return 0;
}