#include <tiffio.h>
#include <iostream>
#include <vector>
#include <cstdint>  // For fixed-width integers

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

void waveletTransform(std::vector<std::uint8_t>& buffer, std::uint32_t width, std::uint32_t height, std::uint16_t samples_per_pixel) {
    // Example: Apply a simple Haar wavelet transform (for illustration).
    for (std::uint32_t row = 0; row < height; row++) {
        for (std::uint32_t col = 0; col < width; col++) {
            buffer[row * width + col] = buffer[row * width + col] / 2;  // Simplified compression
        }
    }
}

void writeTiffImage(const char* filename, const std::vector<std::uint8_t>& buffer, std::uint32_t width, std::uint32_t height, std::uint16_t samples_per_pixel) {
    TIFF* tiff = TIFFOpen(filename, "w");
    if (tiff == nullptr) {
        std::cerr << "Could not open TIFF file for writing!" << std::endl;
        return;
    }

    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, samples_per_pixel);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);  // Assuming 8 bits per channel
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

    for (std::uint32_t row = 0; row < height; row++) {
        TIFFWriteScanline(tiff, const_cast<std::uint8_t*>(&buffer[row * width * samples_per_pixel]), row);
    }

    TIFFClose(tiff);
}

int main() {
    const char* inputFile = "board.tif";
    const char* outputFile = "output_compressed.tif";

    std::uint32_t width, height;
    std::uint16_t samples_per_pixel;
    std::vector<std::uint8_t> buffer;

    // Step 1: Read the TIFF image
    readTiffImage(inputFile, buffer, width, height, samples_per_pixel);

    // Step 2: Apply wavelet compression (simplified)
    waveletTransform(buffer, width, height, samples_per_pixel);

    // Step 3: Write the compressed image to a new TIFF file
    writeTiffImage(outputFile, buffer, width, height, samples_per_pixel);

    std::cout << "Wavelet compression completed!" << std::endl;

    return 0;
}