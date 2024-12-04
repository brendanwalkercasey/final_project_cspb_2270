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

/*
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
*/

// Mexican Hat Wavelet Transform
u32 mexicanHatWavelet(u32 num) {
    return (2.0 / (sqrt(3) * pow(M_PI, 0.25))) * (1 - num * num) * exp(-num * num / 2.0);
}


// Function to apply wavelet compression (simplified)
void waveletTransform(vector<u8>& buffer, u32 width, u32 height) {
    // Example: Apply a simple Haar wavelet transform (for illustration). 
    
    
    // Mexican Hat Wavelet Implimentation
    for (u32 row = 0; row < height; row++){
        for (u32 col = 0; col < width; col++){
            buffer[row * width + col] = mexicanHatWavelet(buffer[row * width + col]);
        }
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
    const char* outputFile = "output_compressed.tif";

    u32 width, height;
    u16 samples_per_pixel;
    vector<u8> buffer;

    // Step 1: Read the TIFF image
    readTiffImage(inputFile, buffer, width, height, samples_per_pixel);

    // Step 2: Convert RGB to Grayscale
    //rgbToGrayscale(buffer, width, height, samples_per_pixel);

    // Step 3: Apply wavelet compression (simplified)
    waveletTransform(buffer, width, height);

    // Step 4: Write the compressed grayscale image to a new TIFF file
    writeTiffImage(outputFile, buffer, width, height, samples_per_pixel);

    cout << "Wavelet compression completed!" << endl;

    return 0;
}