#include <tiffio.h>
#include <iostream>
#include <vector>
#include <cstdint>  // For fixed-width integers

// ***************************************************************************************************************************
/*
Key Changes (from test_libtiff_v02.cpp):
1) rgbToGrayscale function: This function performs the RGB to grayscale conversion using the luminance formula
    - It assumes that the image is in RGB format (3 channels).
    - For each pixel, it calculates the grayscale value and stores it in a new buffer.
2) Writing Grayscale TIFF: When saving the image to a new TIFF file, the pixel format is set to grayscale (PHOTOMETRIC_MINISBLACK and samples_per_pixel = 1).
3) The rest of the code remains mostly the same, except the image is first converted to grayscale before applying wavelet compression.
*/
// ***************************************************************************************************************************


// Function to read the TIFF image into a buffer
void readTiffImage(const char* filename, std::vector<std::uint8_t>& buffer, std::uint32_t& width, std::uint32_t& height, std::uint16_t& samples_per_pixel) {
    TIFF* tiff = TIFFOpen(filename, "r");
    if (tiff == nullptr) {
        std::cerr << "Could not open TIFF file!" << std::endl;
        return;
    }

    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);

    buffer.resize(width * height * samples_per_pixel);
    for (std::uint32_t row = 0; row < height; row++) {
        TIFFReadScanline(tiff, &buffer[row * width * samples_per_pixel], row);
    }

    TIFFClose(tiff);
}

// Function to convert RGB image to grayscale using the luminance formula
void rgbToGrayscale(std::vector<std::uint8_t>& buffer, std::uint32_t width, std::uint32_t height, std::uint16_t samples_per_pixel) {
    // Make a new buffer to store the grayscale image
    std::vector<std::uint8_t> grayscaleBuffer(width * height);

    for (std::uint32_t i = 0; i < height; i++) {
        for (std::uint32_t j = 0; j < width; j++) {
            // Extract the RGB pixel (Assuming 8-bit per channel RGB image)
            std::uint32_t pixel = buffer[i * width * samples_per_pixel + j * samples_per_pixel];

            // Extract RGB components
            std::uint8_t red = (pixel >> 16) & 0xFF;
            std::uint8_t green = (pixel >> 8) & 0xFF;
            std::uint8_t blue = pixel & 0xFF;

            // Compute the grayscale value using the luminance formula
            std::uint8_t grayscale = static_cast<std::uint8_t>(red * 0.299 + green * 0.587 + blue * 0.114);

            // Store the grayscale value in the new buffer (one value per pixel)
            grayscaleBuffer[i * width + j] = grayscale;
        }
    }

    // Copy the grayscale buffer back to the original buffer (if you want to overwrite the original image)
    buffer = std::move(grayscaleBuffer);
}

// Function to apply wavelet compression (simplified)
void waveletTransform(std::vector<std::uint8_t>& buffer, std::uint32_t width, std::uint32_t height) {
    // Example: Apply a simple Haar wavelet transform (for illustration).
    for (std::uint32_t row = 0; row < height; row++) {
        for (std::uint32_t col = 0; col < width; col++) {
            buffer[row * width + col] = buffer[row * width + col] / 2;  // Simplified compression
        }
    }
}

// Function to write the TIFF image from the buffer
void writeTiffImage(const char* filename, const std::vector<std::uint8_t>& buffer, std::uint32_t width, std::uint32_t height, std::uint16_t samples_per_pixel) {
    TIFF* tiff = TIFFOpen(filename, "w");
    if (tiff == nullptr) {
        std::cerr << "Could not open TIFF file for writing!" << std::endl;
        return;
    }

    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);  // Grayscale image will have 1 channel
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);  // 8 bits per pixel for grayscale
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);  // Grayscale, no RGB

    for (std::uint32_t row = 0; row < height; row++) {
        TIFFWriteScanline(tiff, const_cast<std::uint8_t*>(&buffer[row * width]), row);
    }

    TIFFClose(tiff);
}

int main() {
    const char* inputFile = "input.tif";
    const char* outputFile = "output_compressed.tif";

    std::uint32_t width, height;
    std::uint16_t samples_per_pixel;
    std::vector<std::uint8_t> buffer;

    // Step 1: Read the TIFF image
    readTiffImage(inputFile, buffer, width, height, samples_per_pixel);

    // Step 2: Convert RGB to Grayscale
    rgbToGrayscale(buffer, width, height, samples_per_pixel);

    // Step 3: Apply wavelet compression (simplified)
    waveletTransform(buffer, width, height);

    // Step 4: Write the compressed grayscale image to a new TIFF file
    writeTiffImage(outputFile, buffer, width, height, samples_per_pixel);

    std::cout << "Wavelet compression completed!" << std::endl;

    return 0;
}